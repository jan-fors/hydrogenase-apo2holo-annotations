"""
Contains the function create fingerprint, which takes as input:
- structure
- cofactor name/id or point
- radius

returns:
- fingerprint 
"""
from pathlib import Path
from src.utils.constants import (
    FINGERPRINT_RADIUS
)
from Bio.PDB import PDBParser, NeighborSearch
from ase import Atoms
import numpy as np
from ase.io import read
from dscribe.descriptors import SOAP

def create_aminoacid_fingerprint(structure_path : Path, point : tuple, fingerprint_radius : float = FINGERPRINT_RADIUS) -> np.array:
    """"""
    parser = PDBParser()
    structure = parser.get_structure("prot", structure_path)

    atoms = list(structure.get_atoms())
    ns = NeighborSearch(atoms)

    near_atoms = ns.search(point, fingerprint_radius)  

    residues = {a.get_parent() for a in near_atoms}

    F = {
        "ALA": 0,
        "ARG": 0,
        "ASN": 0,
        "ASP": 0,
        "CYS": 0,
        "GLN": 0,
        "GLU": 0,
        "GLY": 0,
        "HIS": 0,
        "ILE": 0,
        "LEU": 0,
        "LYS": 0,
        "MET": 0,
        "PHE": 0,
        "PRO": 0,
        "SER": 0,
        "THR": 0,
        "TRP": 0,
        "TYR": 0,
        "VAL": 0
    }

    for r in residues:
        try:
            F[r.get_resname()] += 1
        except:
            continue
    return np.array(list(F.values()))


def _build_amino_acid_only_structure(structure_path):
    ALLOWED_SPECIES = ["C", "N", "O", "S"]
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", structure_path)

    symbols = []
    positions = []

    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag, resseq, icode = residue.id
                if hetflag != " ":       # " " means standard ATOM record (amino acid)
                    continue             # skip HETATM (cofactors, ions, waters, ligands)
                for atom in residue:
                    element = atom.element.strip()
                    if element not in ALLOWED_SPECIES:
                        continue
                    symbols.append(element)
                    positions.append(tuple(atom.coord))

    return Atoms(symbols=symbols, positions=positions)

def create_physiochemical_radial_angular(structure_path, point, fingerprint_radius=FINGERPRINT_RADIUS):
    
    atoms = _build_amino_acid_only_structure(structure_path)
    ALLOWED_SPECIES = ["C", "N", "O", "S"]
    # keep only atoms whose element is in our declared species list
    keep_mask = [s in ALLOWED_SPECIES for s in atoms.get_chemical_symbols()]
    atoms = atoms[keep_mask]

    soap = SOAP(
        species=ALLOWED_SPECIES,
        r_cut=fingerprint_radius,
        n_max=4,
        l_max=3,
        sigma=0.5,
        periodic=False,
        rbf="gto",
        weighting={"function": "poly", "r0": fingerprint_radius, "c": 1, "m": 3},
    )

    fingerprint = soap.create(atoms, centers=[point])[0]
    
    return np.array(fingerprint)

#########

# --- Per-residue physicochemical property lookup tables -------------------

# Kyte-Doolittle hydrophobicity scale (higher = more hydrophobic)
# Source: Kyte J., Doolittle R.F. (1982). "A simple method for displaying
# the hydropathic character of a protein." J. Mol. Biol. 157(1):105-132.
# https://doi.org/10.1016/0022-2836(82)90515-0
HYDROPHOBICITY = {
    "ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5,
    "GLN": -3.5, "GLU": -3.5, "GLY": -0.4, "HIS": -3.2, "ILE": 4.5,
    "LEU": 3.8, "LYS": -3.9, "MET": 1.9, "PHE": 2.8, "PRO": -1.6,
    "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2,
}

# Formal side-chain charge at physiological pH
CHARGE = {
    "ASP": -1.0, "GLU": -1.0,
    "LYS": 1.0, "ARG": 1.0, "HIS": 0.1,
}

# Residues with aromatic side chains
AROMATIC = {"PHE", "TYR", "TRP", "HIS"}

# Residues capable of side-chain hydrogen bond donation / acceptance
HBOND_DONOR = {"SER", "THR", "TYR", "ASN", "GLN", "LYS", "ARG", "HIS", "TRP", "CYS"}
HBOND_ACCEPTOR = {"SER", "THR", "TYR", "ASN", "GLN", "ASP", "GLU", "HIS"}


def _residue_properties(resname):
    """Return a fixed-order property vector for one residue type."""
    return np.array([
        HYDROPHOBICITY.get(resname, 0.0),
        CHARGE.get(resname, 0.0),
        1.0 if resname in AROMATIC else 0.0,
        1.0 if resname in HBOND_DONOR else 0.0,
        1.0 if resname in HBOND_ACCEPTOR else 0.0,
    ])


N_PROPERTIES = 6  # hydrophobicity, charge, aromatic, donor, acceptor, sulfur


def create_feature_fingerprint(structure_path, point, fingerprint_radius=7.5, n_shells=6):
    """
    Build a FEATURE-style fingerprint around a 3D point.

    Parameters
    ----------
    structure_path : str or Path
        Path to the PDB file.
    point : tuple(float, float, float)
        Center of the microenvironment (e.g. active site / pocket center).
    fingerprint_radius : float
        Outer radius of the microenvironment, in Angstrom.
    n_shells : int
        Number of concentric shells to divide the radius into.

    Returns
    -------
    np.ndarray of shape (n_shells * N_PROPERTIES,)
        Per-shell summed physicochemical properties. Shell 0 is innermost.
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", structure_path)

    point = np.array(point, dtype=float)
    shell_edges = np.linspace(0, fingerprint_radius, n_shells + 1)
    fingerprint = np.zeros(n_shells * N_PROPERTIES)

    for model in structure:
        for chain in model:
            for residue in chain:
                hetflag, resseq, icode = residue.id
                if hetflag != " ":
                    continue  # skip HETATMs (cofactors, ions, waters, ligands)

                resname = residue.resname
                residue_props = _residue_properties(resname)  # length 5

                for atom in residue:
                    dist = np.linalg.norm(atom.coord - point)
                    if dist >= fingerprint_radius:
                        continue

                    shell_idx = np.searchsorted(shell_edges, dist, side="right") - 1
                    shell_idx = min(shell_idx, n_shells - 1)

                    # Sulfur is atom-specific (e.g. Cys S-gamma, Met S-delta),
                    # not shared by every atom in the residue -- checked per atom.
                    is_sulfur = 1.0 if atom.element.strip().upper() == "S" else 0.0
                    atom_props = np.append(residue_props, is_sulfur)  # length 6

                    start = shell_idx * N_PROPERTIES
                    fingerprint[start:start + N_PROPERTIES] += atom_props

    return fingerprint

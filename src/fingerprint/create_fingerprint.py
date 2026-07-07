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
    FINGERPRINT_RADIUS,
    SOAP_N_MAX,
    SOAP_L_MAX,
    SOAP_ALLOWED_SPECIES,
    SOAP_RBF,
    SOAP_SIGMA,
)
from Bio.PDB import PDBParser, NeighborSearch
from ase import Atoms
import numpy as np
from ase.io import read
from dscribe.descriptors import SOAP


def create_aminoacid_fingerprint(
    structure_path: Path, point: tuple, fingerprint_radius: float = FINGERPRINT_RADIUS
) -> np.array:
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
        "VAL": 0,
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
                if hetflag != " ":  # " " means standard ATOM record (amino acid)
                    continue  # skip HETATM (cofactors, ions, waters, ligands)
                for atom in residue:
                    element = atom.element.strip()
                    if element not in ALLOWED_SPECIES:
                        continue
                    symbols.append(element)
                    positions.append(tuple(atom.coord))

    return Atoms(symbols=symbols, positions=positions)


def create_physiochemical_radial_angular(
    structure_path, point, fingerprint_radius=FINGERPRINT_RADIUS
):

    atoms = _build_amino_acid_only_structure(structure_path)
    keep_mask = [s in SOAP_ALLOWED_SPECIES for s in atoms.get_chemical_symbols()]
    atoms = atoms[keep_mask]

    soap = SOAP(
        species=SOAP_ALLOWED_SPECIES,
        r_cut=fingerprint_radius,
        n_max=SOAP_N_MAX,
        l_max=SOAP_L_MAX,
        sigma=SOAP_SIGMA,
        periodic=False,
        rbf=SOAP_RBF,
        weighting={"function": "poly", "r0": fingerprint_radius, "c": 1, "m": 3},
    )

    fingerprint = soap.create(atoms, centers=[point])[0]

    return np.array(fingerprint)


def create_combined_fingerprint(
    structure_path, point, fingerprint_radius=FINGERPRINT_RADIUS
):
    """ """
    aa_fp = create_aminoacid_fingerprint(structure_path, point, fingerprint_radius)
    pra_fp = create_physiochemical_radial_angular(
        structure_path, point, fingerprint_radius
    )

    return np.concatenate([aa_fp, pra_fp])

"""
Contains the function create fingerprint, which takes as input:
- structure
- cofactor name/id or point
- radius

returns:
- fingerprint
"""

from pathlib import Path
from src.pdb.pdb_handler import get_atom_positions_for_prca, get_structure
import numpy as np
from ase.io import read
from dscribe.descriptors import SOAP
from Bio.PDB.NeighborSearch import NeighborSearch
from src.utils.constants import STANDARD_AMINO_ACIDS
from typing import List
import tqdm


def count_aminoacids(
        structure_path, 
        point, 
        f_radius,
        params
) -> np.array:
    """
    """
    structure = get_structure(structure_path)
    atoms = list(structure.get_atoms())
    ns = NeighborSearch(atoms)
    near_atoms = ns.search(point, f_radius)
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

def count_aminoacids_per_dist(
        structure_path, 
        point, 
        f_radius,
        params):
    """
    """
    point = np.asarray(point, dtype=float)
 
    shell_width = params["shell_width"]
    n_shells = int(np.ceil(f_radius / shell_width))
    fingerprint = np.zeros(20 * n_shells, dtype=int)
    elem_index = {e: i for i, e in enumerate(STANDARD_AMINO_ACIDS)}
 
    structure = get_structure(structure_path)
 
    to_remove = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":  
                    to_remove.append((chain, residue.id))
 
    for chain, res_id in to_remove:
        chain.detach_child(res_id)
 
    residues = []
    for model in structure:
        for chain in model:
            for residue in chain:
                residues.append((model.id, chain.id, residue.id, residue))
    residues.sort(key=lambda r: (r[0], r[1], r[2]))
 
    for _, _, _, residue in residues:
        resname = residue.get_resname().strip().upper()
        if resname not in elem_index:
            continue
 
        atoms = residue.get_unpacked_list()
        if not atoms:
            continue
        dist = min(np.linalg.norm(atom.coord - point) for atom in atoms)
 
        if dist >= f_radius:
            continue
 
        shell = int(dist // shell_width)  # 0 for [0, shell_width), 1 for next, ...
        idx = shell * len(STANDARD_AMINO_ACIDS) + elem_index[resname]
        fingerprint[idx] += 1
 
    return fingerprint


def physiochemical_radial_angular(
    structure_path, point, fingerprint_radius, params
):
    """
    """
    allowed_species = params["allowed_species"]
    n_max = params["n_max"]
    l_max = params["l_max"]
    sigma = params["sigma"]
    rbf = params["rbf"]

    atoms = get_atom_positions_for_prca(structure_path, allowed_species)
    keep_mask = [s in allowed_species for s in atoms.get_chemical_symbols()]
    atoms = atoms[keep_mask]

    atoms = atoms[keep_mask]
    
    soap = SOAP(
        species=allowed_species,
        r_cut=fingerprint_radius,
        n_max=n_max,
        l_max=l_max,
        sigma=sigma,
        periodic=False,
        rbf=rbf,
        weighting={"function": "poly", "r0": fingerprint_radius, "c": 1, "m": 3},
    )

    fingerprint = soap.create(atoms, centers=[point])[0]

    return np.array(fingerprint)

def count_atoms_per_dist(
    structure_path: Path, point: tuple, fingerprint_radius: float, params
) -> np.array:
    """
    """
    point = np.asarray(point, dtype=float)

    shell_width = params["shell_width"]
    allowed_species = params["allowed_species"]

    n_shells = int(np.ceil(fingerprint_radius / shell_width))
    fingerprint = np.zeros(len(allowed_species) * n_shells, dtype=int)
    elem_index = {e: i for i, e in enumerate(allowed_species)}

    structure = get_structure(structure_path)

    # collect hetero residues first, then detach (don't mutate while iterating)
    to_remove = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != " ":        # non-blank hetflag = HETATM (incl. water)
                    to_remove.append((chain, residue.id))

    for chain, res_id in to_remove:
        chain.detach_child(res_id)

    for atom in structure.get_atoms():
        dist = np.linalg.norm(atom.coord - point)
        if dist >= fingerprint_radius:
            continue

        element = atom.element.strip().capitalize()
        if element not in elem_index:
            continue

        shell = int(dist // shell_width)  # 0 for [0,1), 1 for [1,2), ...
        idx = shell * len(allowed_species) + elem_index[element]
        fingerprint[idx] += 1

    return np.array(fingerprint)

def count_atoms(structure_path: Path, point: tuple, fingerprint_radius: float, params):
    """
    """
    structure = get_structure(structure_path)
    atoms = list(structure.get_atoms())
    ns = NeighborSearch(atoms)
    near_atoms = ns.search(point, fingerprint_radius)
    
    F = {}

    for atom_type in params["allowed_species"]:
        F[atom_type] = 0

    for a in near_atoms:
        try:
            F[a.get_name()()] += 1
        except:
            continue
    return np.array(list(F.values()))


FINGERPRINT_REGISTRY = {
    "aa_count": count_aminoacids,
    "aa_dist": count_aminoacids_per_dist,
    "ac_count": count_atoms,
    "ac_dist": count_atoms_per_dist,
    "pcra": physiochemical_radial_angular
}

def create_fingerprint(structure_path : Path, point : tuple[float,float,float], fingerprint_radius : float, fingerprint_types : list):
    """
    """
    # sort fingerprint_types
    fingerprint_types_sorted = sorted(fingerprint_types, key=lambda x: x["name"])

    fps = []

    for fp in fingerprint_types_sorted:
        fp_type = fp["name"]
        fp_fn = FINGERPRINT_REGISTRY[fp_type]
        params = fp["params"]
        
        fps.append(fp_fn(structure_path, point, fingerprint_radius, params))

    return np.concatenate(fps)

def create_fingerprints(input_path: Path,
        coordinates: List[tuple[float, float, float]],
        f_radius: float,
        fingerprint_types : list
    ):
        """"""
        fingerprints = []
        for point in tqdm(coordinates):
            # create fingerprint
            F = create_fingerprint(input_path, point, f_radius, fingerprint_types)
            fingerprints.append(F)
    
        return fingerprints
"""
"""
# imports
from apo2holo.pdb.pdb_handler import get_structure
import numpy as np
from Bio.PDB.NeighborSearch import NeighborSearch

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
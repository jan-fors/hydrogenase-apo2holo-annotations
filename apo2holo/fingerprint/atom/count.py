"""
"""
# imports
from apo2holo.pdb.pdb_handler import get_structure
from pathlib import Path
from Bio.PDB.NeighborSearch import NeighborSearch
import numpy as np

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
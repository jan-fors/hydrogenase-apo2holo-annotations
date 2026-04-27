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
import numpy as np

def create_fingerprint(structure_path : Path, point : tuple) -> np.array:
    """"""
    parser = PDBParser()
    structure = parser.get_structure("prot", structure_path)

    atoms = list(structure.get_atoms())
    ns = NeighborSearch(atoms)

    near_atoms = ns.search(point, FINGERPRINT_RADIUS)  

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
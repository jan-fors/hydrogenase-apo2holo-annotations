"""
Contains the function create fingerprint, which takes as input:
- structure
- cofactor name/id or point
- radius

returns:
- fingerprint
"""

from pathlib import Path

import numpy as np

from apo2holo.utils.constants import STANDARD_AMINO_ACIDS
from typing import List
from tqdm import tqdm
from apo2holo.fingerprint.aminoacid.count import count_aminoacids
from apo2holo.fingerprint.aminoacid.shell_count import count_aminoacids_per_dist
from apo2holo.fingerprint.atom.shell_count import count_atoms_per_dist
from apo2holo.fingerprint.atom.count import count_atoms
from apo2holo.fingerprint.graph.wl_vec import create_wl_hash


FINGERPRINT_REGISTRY = {
    "aa_count": count_aminoacids,
    "aa_dist": count_aminoacids_per_dist,
    "ac_count": count_atoms,
    "ac_dist": count_atoms_per_dist,
    "wl_hash": create_wl_hash
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
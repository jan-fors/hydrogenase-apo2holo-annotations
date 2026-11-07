"""
"""
# imports
from apo2holo.utils.constants import STANDARD_AMINO_ACIDS
from apo2holo.pdb.pdb_handler import get_structure
import numpy as np

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
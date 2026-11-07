"""
"""
# imports
from pathlib import Path
import numpy as np
from apo2holo.pdb.pdb_handler import get_structure

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
from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch, Selection
import numpy as np

def identify_cysteines(structure_path : str, coords : list, radius : float = 6.0):
        """
    Identify all cysteine residues within a given radius of one or more coordinates.

    Parameters
    ----------
    structure_path : str
        Path to a PDB or mmCIF structure file.
    coords : list
        List of 3D coordinates, e.g. [(x, y, z), (x, y, z)].
    radius : float, default=4.0
        Search radius in Angstrom.

    Returns
    -------
    list
        List of tuples: [(chain_id, residue_number), ...]
    """
        parser = PDBParser(QUIET=True)

        structure = parser.get_structure("structure", structure_path)

        atoms = Selection.unfold_entities(structure, "A")
        ns = NeighborSearch(atoms)

        # wichtig: als numpy array
        center = np.array(coords, dtype=float)

        found_cys = set()

        nearby_atoms = ns.search(center, radius, level="A")

        for atom in nearby_atoms:
            residue = atom.get_parent()
            if residue.get_resname().strip() == "CYS":
                chain_id = residue.get_parent().id
                residue_number = residue.id[1]
                found_cys.add((chain_id, residue_number))

        return sorted(found_cys, key=lambda x: (x[0], x[1]))
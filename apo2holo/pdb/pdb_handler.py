from ase import Atoms
from Bio.PDB import PDBParser, NeighborSearch
from pathlib import Path
from apo2holo.io.writers.printl import printl
from typing import List
import os
from Bio.PDB import PDBParser, MMCIFParser, NeighborSearch, Selection
from Bio import SeqIO
import subprocess
import uuid
import numpy as np

def extract_sequences(structure_path) -> dict:
    """
    """
    res = {}
    for record in SeqIO.parse(structure_path, "pdb-atom"):
        res[record.annotations["chain"]] = str(record.seq)

    return res

def get_atom_positions_for_prca(structure_path, allowed_species):
    """
    """
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
                    if element not in allowed_species:
                        continue
                    symbols.append(element)
                    positions.append(tuple(atom.coord))

    return Atoms(symbols=symbols, positions=positions)

def get_structure(structure_path : Path):
    """
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", str(structure_path))

    return structure

def extract_protein_chains_from_file(input_path: Path, out: Path) -> List[Path]:
    """ """
    # extract chains
    printl("Extract chains from protein file")
    chains = get_chains(input_path)
    chain_paths = []
    for chain in chains:
        chain_paths.append(extract_chain(input_path, out, chain))

    return chain_paths

def get_chains(input_pdb: str) -> List[str]:
    """ """
    chains = set()
    with open(input_pdb, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 22:
                c = line[21].strip()
                if c:
                    chains.add(c)
    return sorted(chains)

def extract_chain(input_pdb: str, output_dir: str, chain: str):
    """ """
    structure_name = os.path.basename(input_pdb).split(".")[0] + f"_{chain}.pdb"
    outpath = Path(os.path.join(output_dir, structure_name))
    cmd = ["pdb_selchain", f"-{chain}", str(input_pdb)]

    with outpath.open("w") as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True, check=True)

    return outpath

def extract_cofactors(pdb_path: str) -> dict[str, list[tuple[str, float, float, float]]]:
    """
    Returns:
        {
            "RESNAME": [(element, x, y, z), ...],
            ...
        }
    """

    het_groups = {}

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("HETATM"):

                resname = line[17:20].strip()
                element = line[76:78].strip()

                # Fallback: falls Element-Spalte leer ist
                if not element:
                    element = line[12:16].strip()[0]

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                if resname not in het_groups:
                    het_groups[resname] = []

                het_groups[resname].append((element, x, y, z))

    return het_groups

def identify_cysteines(structure_path : str, coords : list, radius : float = 5.0):
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

def get_protein_boundaries(structure_path: Path) -> tuple:
    """
    """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("s", structure_path)

    coords = np.array([atom.coord for atom in structure.get_atoms()])
    mins = coords.min(axis=0)
    maxs = coords.max(axis=0)

    x_min, x_max, y_min, y_max, z_min, z_max = (
        mins[0],
        maxs[0],
        mins[1],
        maxs[1],
        mins[2],
        maxs[2],
    )
    return x_min, x_max, y_min, y_max, z_min, z_max


def extract_hetatm_residues(pdb_file, exclude_water=True):
    """ """
    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("struct", pdb_file)

    results = {}
    for model in structure:
        for chain in model:
            for residue in chain:

                hetflag, resseq, icode = residue.id

                if not hetflag.startswith("H_"):
                    continue

                if exclude_water and residue.resname in {"HOH", "WAT"}:
                    continue

                atoms = []
                for atom in residue:
                    element = atom.element.strip()
                    if (
                        element == "X"
                        and atom.get_name().strip().upper().startswith("FE")
                    ):
                        element = "FE"
                    x, y, z = atom.coord
                    atoms.append((element, float(x), float(y), float(z)))

                unique_identifier = uuid.uuid4()

                results[unique_identifier] = {
                    "res_name": residue.resname,
                    "res_id": resseq,
                    "chain": chain.id,
                    "atoms": atoms,
                }

    return dict(results)
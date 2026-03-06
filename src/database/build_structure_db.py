"""
Takes a folder which has PDB files in it and creates a database from it.

To achieve this:
1. extract sequences and create seq file
2. create a specific folder structure to search fro homology search
3. create fingerprint db
"""

import argparse
import subprocess
import os
from typing import List
from src.io.printl import printl
from pathlib import Path
from src.utils.constants import COFACTOR_BLACKLIST, STRUCTURE_DB
from src.utils.get_chains import get_chains
from src.utils.extract_chain import extract_chain


def _check_input_dir(input_dir: str):
    pass


def _check_output_dir(output_dir: str):
    """ """
    if not os.path.exists(output_dir):
        printl(f"{output_dir} does not exist. Creating it ...")
        os.makedirs(output_dir, exist_ok=True)


def build_structure_db(input_dir: str, output_dir: str):
    """ """
    _check_input_dir(input_dir)
    _check_output_dir(output_dir)

    for file in os.listdir(input_dir):
        """ """
        file_path = os.path.join(input_dir, file)
        if not os.path.exists(file_path):
            printl(f"{file_path} does not exist.")
            continue

        # get chains
        chains = get_chains(file_path)

        # split each structure into subunits
        for chain in chains:
            result_structure = extract_chain(file_path, output_dir, chain)
            _apply_blacklist(result_structure)

    _create_foldseek_db(output_dir)

def _apply_blacklist(structure_path: str):
    """ """
    pdb_path = Path(structure_path)
    tmp_path = pdb_path.with_suffix(structure_path.suffix + ".tmp")

    removed = 0

    with pdb_path.open("r") as fin, tmp_path.open("w") as fout:
        for line in fin:
            record = line[0:6].strip()
            if record == "HETATM":
                resname = line[17:20].strip().upper()
                if resname in COFACTOR_BLACKLIST:
                    removed += 1
                    continue
            fout.write(line)

    # Atomar ersetzen (sehr wichtig!)
    os.replace(tmp_path, pdb_path)

    printl(f"[remove_blacklisted_hetatm_inplace] removed_hetatm_lines={removed}")


def _create_foldseek_db(structure_dir_path: str):
    """"""
    # create db
    cmd = ["foldseek", "createdb", str(structure_dir_path), STRUCTURE_DB]
    subprocess.run(cmd, stderr=subprocess.PIPE, text=True, check=True)

    # create index
    cmd = ["foldseek", "createindex", STRUCTURE_DB, "tmp"]
    subprocess.run(cmd, stderr=subprocess.PIPE, text=True, check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", help="Path to the structure folder.")
    parser.add_argument("output_dir", help="Path to the structure files.")
    args = parser.parse_args()
    build_structure_db(args.input_dir, args.output_dir)

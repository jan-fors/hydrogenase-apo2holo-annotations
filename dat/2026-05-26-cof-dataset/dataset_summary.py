"""
- What cofactors are in the dataset?
- In what abundance?
"""
import os
from pathlib import Path
import argparse
from collections import Counter
from Bio import PDB

WHITELIST = ["F3S", "F4S", "SF3", "SF4"]

def summary(directory : Path):
    """
    iterate over each file in the provided directory and count the cofactors
    """
    hetatm_counter = Counter()  
    counter = 0

    parser = PDB.PDBParser(QUIET=True)
    for file in os.listdir(directory):
        structure = parser.get_structure("my_protein", directory / Path(file))

        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.get_id()[0] != " ":
                        hetatm_counter[residue.get_resname().strip()] += 1
            break
        counter += 1

    print("Structures", str(counter))

    print("="*10, "FeS Cluster", "="*10)
    for resname, count in hetatm_counter.most_common():
        if resname in WHITELIST:
            print(f"{resname}: {count}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("dir", type=Path)

    args = parser.parse_args()

    summary(args.dir)
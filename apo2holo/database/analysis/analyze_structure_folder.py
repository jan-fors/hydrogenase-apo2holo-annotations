"""
python -m apo2holo.database.analysis.analyze_structure_folder dat/2026-05-26-cof-dataset/seq_repr_v3
"""

import argparse
from pathlib import Path
import os
from apo2holo.pdb.pdb_handler import extract_cofactors
from collections import Counter
from pprint import pprint
from apo2holo.utils.constants import COFACTOR_WHITELIST

def analyze_structure_folder(directory : Path):
    """
    """
    n_stuctures = len(os.listdir(directory))
    print("n_structures", n_stuctures)

    print("\nCofactors in structures")
    counter = Counter([])
    for file in os.listdir(directory):
        cofactors = extract_cofactors(directory/file)
        counter.update(list(cofactors.keys()))  

        #print(cofactors.keys())


    counter = dict(counter)
    
    negatives = 0
    for key in counter.keys():
        if key in COFACTOR_WHITELIST:
            print(f"{key} : {counter[key]}")
        else:
            negatives += 1
    print(f"Negatives: {negatives}")

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()

    analyze_structure_folder(args.directory)

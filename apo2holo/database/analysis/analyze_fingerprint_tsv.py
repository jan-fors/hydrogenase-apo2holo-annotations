"""
python -m apo2holo.database.analysis.analyze_fingerprint_tsv benchmarking/complete_benchmark/2026_09_01_16_16_experiment_04/0/db/fingerprintDB/dd_experiment_04_fpr7.0.tsv
"""

import argparse
from pathlib import Path
import pandas as pd
from collections import Counter
from pprint import pprint
from apo2holo.utils.constants import COFACTOR_WHITELIST

def analyze_fingerprint_tsv(fingerprint_tsv : Path):
    """
    """
    df = pd.read_csv(fingerprint_tsv, sep="\t")

    n_entries = df.shape[0]
    print(f"n_entries: {n_entries}")

    
    try:
        print("\nres_name distribution:")
        res_name_list = df["res_name"].to_list()
        counter = dict(Counter(res_name_list))

        negatives = 0
        for key in counter.keys():
            if key in COFACTOR_WHITELIST:
                print(f"{key} : {counter[key]}")
            else:
                negatives += 1
        print(f"Negatives: {negatives}")
    except:
        print("No Resname")

    print("\nformula distribution")
    formula_list = df["formula"].to_list()
    counter = dict(Counter(formula_list))

    pprint(counter)

    print("\ntype distribution:")
    type_list = df["type"].to_list()
    counter = dict(Counter(type_list))

    pprint(counter)
   
    #print(df.columns)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("fingerprint_tsv",type=Path)
    args = parser.parse_args()

    analyze_fingerprint_tsv(args.fingerprint_tsv)
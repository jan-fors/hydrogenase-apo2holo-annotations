"""
Reads the result of a foldseek run and copies the representatives into a new folder
"""
import os
from pathlib import Path
import argparse
import pandas as pd
from typing import List
import shutil

def extract_unique_structures(foldseek_result_tsv : Path, src_dirs : List[Path], dest_dir : Path):
    """
    """
    # read foldseek_result_tsv
    df = pd.read_csv(foldseek_result_tsv, sep="\t", header=None)
    df.columns = ["repr", "samp"]
    
    # get representatives
    unique_repr = list(df["repr"].unique())
    unique_repr = [x.split("_")[0] for x in unique_repr]
    
    for i in unique_repr:
        dest_path = dest_dir / Path(i + ".pdb")
        
        if os.path.exists(dest_path):
            continue

        for j in src_dirs:
            src_path = j / Path(i + ".pdb")

            if os.path.exists(src_path):
                shutil.copy2(src_path, dest_path)

    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("foldseek_res_tsv", type=Path, help="Path to the foldseek result tsv file")
    parser.add_argument("--src_dirs", nargs="+", help="Path to the src directories", required=True)
    parser.add_argument("--dest_dir", type=Path, default="out", help="Destiation")

    args = parser.parse_args()

    # validate input
    if not os.path.exists(args.foldseek_res_tsv):
        raise ValueError(f"{args.foldseek_res_tsv} does not exist")
    
    for i in args.src_dirs:
        if not os.path.exists(i):
            raise ValueError(f"{i} does not exist")
        
    if not os.path.exists(args.dest_dir):
        os.makedirs(args.dest_dir, exist_ok=True)

    extract_unique_structures(args.foldseek_res_tsv, args.src_dirs, args.dest_dir)
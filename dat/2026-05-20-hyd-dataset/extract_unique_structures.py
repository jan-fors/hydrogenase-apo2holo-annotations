"""
Read result tsv and extract all representatives into a defined folder
"""
import os
import argparse
from pathlib import Path
import pandas as pd
import shutil

def extract_unique_structures(foldseek_result_tsv : Path, source_folder : Path, dest_folder : Path):
    """
    
    """
    # read
    df = pd.read_csv(foldseek_result_tsv, sep="\t", header=None)
    df.columns = ["representatives", "samples"]

    # get unique repr
    representatives = df["representatives"].to_list()
    representatives = list(set(representatives))
    
    for r in representatives:
        src = source_folder / Path(r + ".pdb")
        dest = dest_folder / Path(r.split("_")[0] + ".pdb")
        
        shutil.copy2(src,dest)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    
    parser.add_argument("foldseek_result_tsv", help="Path to the result tsv table.")
    parser.add_argument("source_folder", help="Path to the source folder")
    parser.add_argument("dest_folder", help="Path to the dest folder")

    args = parser.parse_args()


    if not os.path.exists(args.foldseek_result_tsv):
        raise ValueError("Provide valid foldseek tsv path")
    

    if not os.path.exists(args.source_folder):
        raise ValueError("Provide valid source_folder path")
    

    if not os.path.exists(args.dest_folder):
        raise ValueError("Provide valid dest_folder path")

    extract_unique_structures(args.foldseek_result_tsv, args.source_folder, args.dest_folder)

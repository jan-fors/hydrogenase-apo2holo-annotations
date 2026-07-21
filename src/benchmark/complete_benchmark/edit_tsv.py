import argparse
import pandas as pd
from pathlib import Path
import os

REMOVE_TYPE = None
FP = "atom_count"
AUG_DIST = None

def edit_tsv(dir: Path, data_path_og : Path, out_name : str):
    """
    """
    for subset in os.listdir(dir):
        data_path = dir / Path(subset) / data_path_og
        # read data
        df = pd.read_csv(data_path, sep="\t")
        print(df.shape)
        print(df.head())

        if REMOVE_TYPE != None:
            # remove type from tsv
            rows = df[df["type"] == REMOVE_TYPE].index
            df.drop(rows, inplace=True)
            
            print(f"Removed type {REMOVE_TYPE}")
            print(df.shape)
        
        if FP != None:
            if FP == "atom_count":
                # remove the first 20 columns
                df = df.drop(df.columns[range(0,20)], axis=1)
                print("Removed Columns")
                print(df.shape)
                
            if FP == "aminoacid":
                # remove the 20 to 43 columns
                df = df.drop(df.columns[range(20,43)], axis=1)
                print("Removed Columns")
                print(df.shape)

        if AUG_DIST != None:
            rows = df[df["aug_dist"] > AUG_DIST].index
            df = df.drop(rows, inplace=True)
            
            print(f"Removed aug dist higher than {AUG_DIST}")
            print(df.shape)

        df.to_csv(data_path.parent/ Path(out_name+".tsv"), sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=Path)
    parser.add_argument("data", type=Path)
    parser.add_argument("out_name", type=str)
    args = parser.parse_args()

    edit_tsv(args.dir, args.data, args.out_name)

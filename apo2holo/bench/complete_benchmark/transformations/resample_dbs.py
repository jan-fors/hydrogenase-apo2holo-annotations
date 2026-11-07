"""
python -m src.benchmark.complete_benchmark.resample_dbs benchmarking/complete_benchmark/subsets
"""
import argparse
import os
import pandas as pd
from pathlib import Path
from apo2holo.database.utils.random_sampler import random_oversample

def resample_dbs(directory_path: Path, db_name : str, undersampling : bool, top : int) -> str:
    """
    Method that resamples databases specified by name.
    Iterates over all subsets and applies resampling rules.
    Return name of updated database
    """
    for subset in os.listdir(directory_path):
        # define db path
        db_path = directory_path/Path(subset)/Path('db')/Path('fingerprintDB')/Path(db_name)
        
        # read db
        if str(db_path).endswith(".gz"):
            df = pd.read_csv(db_path, sep="\t", compression="gzip")
        else:
            df = pd.read_csv(db_path, sep="\t")

        # extract class column -> formula
        Y = df["formula"]
        try:
            # drop all non-fingeprint columns
            X = df.drop(columns=['formula', 'id', 'smiles', 'res_name'])
        except Exception as e:
            print(e)
            try:
                 X = df.drop(columns=['formula'])
            except Exception as e:
                print(e)

        # perform oversampling
        new_X, new_Y = random_oversample(X, Y, target_count=top, undersample=undersampling)
        new_X = pd.DataFrame(new_X, columns=X.columns)   # X = the frame before resampling
        new_Y = pd.Series(new_Y, name="formula")
        new_df = pd.concat([new_Y, new_X], axis=1)

        # define output path
        output_file_name = f"r{top}_U{undersampling}_"+db_name
        out_path = directory_path/Path(subset)/Path('db')/Path('fingerprintDB')/Path(output_file_name)

        # save
        if str(out_path).endswith(".gz"):
            new_df.to_csv(out_path, index=False, sep="\t", compression="gzip")
        else:
            new_df.to_csv(out_path, index=False, sep="\t")

    return output_file_name

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=Path)
    parser.add_argument("db_name", type=str)

    args = parser.parse_args()

    resample_dbs(args.dir, args.db_name)
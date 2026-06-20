"""
python -m src.benchmark.complete_benchmark.resample_dbs benchmarking/complete_benchmark/subsets
"""
import argparse
import os
import pandas as pd
from pathlib import Path
from src.database.random_sampler import random_oversample


def resample_dbs(directory_path: Path, db_name : str):
    """
    """
    for subset in os.listdir(directory_path):
        db_path = directory_path/Path(subset)/Path('db')/Path('fingerprintDB')/Path(db_name)
        
        # read db
        df = pd.read_csv(db_path, sep="\t")
        
        Y = df["formula"]
        try:
            X = df.drop(columns=['formula', 'id', 'smiles', 'res_name'])
        except Exception as e:
            print(e)
            try:
                 X = df.drop(columns=['formula'])
            except Exception as e:
                print(e)

        new_X, new_Y = random_oversample(X, Y, target_count=1000, undersample=True)
        new_X = pd.DataFrame(new_X, columns=X.columns)   # X = the frame before resampling
        new_Y = pd.Series(new_Y, name="formula")
        out = pd.concat([new_Y, new_X], axis=1)

        out_path = directory_path/Path(subset)/Path('db')/Path('fingerprintDB')/Path("resampled_"+db_name)
        out.to_csv(out_path, index=False, sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=Path)
    parser.add_argument("db_name", type=str)

    args = parser.parse_args()

    resample_dbs(args.dir, args.db_name)
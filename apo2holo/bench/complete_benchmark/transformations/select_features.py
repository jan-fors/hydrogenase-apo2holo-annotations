import argparse
import os
import pandas as pd
from pathlib import Path
from apo2holo.database.utils.feature_selector import feature_selector


def select(directory_path: Path, db_name : str, keep:list):
    """
    """
    for subset in os.listdir(directory_path):
        # define db path
        db_path = directory_path/Path(subset)/Path('db')/Path('fingerprintDB')/Path(db_name)
        
        # read db
        if str(db_path).endswith(".gz"):
            df = pd.read_csv(db_path, sep="\t", compression="gzip")
        else:
            df = pd.read_csv(db_path, sep="\t")

        # select features
        new_df = feature_selector(df, keep)

        # define output path
        output_file_name = "sf_"+db_name
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
    parser.add_argument("--keep", nargs="+", help="one or more strings")

    args = parser.parse_args()

    select(args.dir, args.db_name, args.keep)
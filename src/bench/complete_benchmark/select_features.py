import argparse
import os
import pandas as pd
from pathlib import Path
from src.database.utils.feature_selector import feature_selector


def select(directory_path: Path, db_name : str, keep:list):
    """
    """
    for subset in os.listdir(directory_path):
        db_path = directory_path/Path(subset)/Path('db')/Path('fingerprintDB')/Path(db_name)
        
        # read db
        df = pd.read_csv(db_path, sep="\t")

        out = feature_selector(df, keep)
        
        out_path = directory_path/Path(subset)/Path('db')/Path('fingerprintDB')/Path("sf_"+db_name)
        out.to_csv(out_path, index=False, sep="\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=Path)
    parser.add_argument("db_name", type=str)
    parser.add_argument("--keep", nargs="+", help="one or more strings")

    args = parser.parse_args()

    select(args.dir, args.db_name, args.keep)
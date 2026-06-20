"""
python -m src.benchmark.complete_benchmark.deduplicate_fingeprint_db benchmarking/complete_benchmark/subsets db_R5.0.tsv
"""
import argparse
import os
import pandas as pd
from pathlib import Path
from src.database.deduplicator import deduplicate

def deduplicate_fingerprint_db(directory_path : Path, db_name : str):
    for subset in os.listdir(directory_path):
        db_path = directory_path/Path(subset)/Path('db')/Path('fingerprintDB')/Path(db_name)
        
        # read db
        df = pd.read_csv(db_path, sep="\t")

        new_df = deduplicate(df)

        out_path = directory_path/Path(subset)/Path('db')/Path('fingerprintDB')/Path("dd_"+db_name)
        new_df.to_csv(out_path, index = False, sep = "\t")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=Path)
    parser.add_argument("db_name", type=str)

    args = parser.parse_args()

    deduplicate_fingerprint_db(args.dir, args.db_name)
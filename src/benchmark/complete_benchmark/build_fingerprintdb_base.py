"""
creates the fingerprint.tsv in the correct folder
"""
import os
from pathlib import Path
import argparse
from src.database.FingerprintDB import FingerprintDB

def build_fingerprint_base(directory : Path, f_radius : float):
    """
    """
    # iterate over each subset
    for subset in os.listdir(directory):
        print("Build fingerprint table for", subset)
        fDB = FingerprintDB()

        subset_path = directory / Path(str(subset))

        train_dir = subset_path / Path('train')
        db_path = subset_path / Path('db')

        fingerprint_db_path = db_path / Path('fingerprintDB')
        os.makedirs(fingerprint_db_path, exist_ok=True)

        out_file = fingerprint_db_path / Path("db_R" + str(f_radius) + ".tsv")

        fDB.build_db(train_dir, f_radius, False)
        fDB.save_fingerprint_tsv(out_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("dir", type=Path)
    parser.add_argument("--f_radius", type=float, default=5.0)
    args = parser.parse_args()

    build_fingerprint_base(args.dir, args.f_radius)
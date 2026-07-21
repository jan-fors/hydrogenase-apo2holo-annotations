"""
creates the fingerprint.tsv in the correct folder
"""

import os
from pathlib import Path
import argparse
from src.database.FingerprintDB import FingerprintDB


def build_fingerprint_base(directory: Path, f_radius: float, jobs: int, name : str):
    """
    input_directory: Path,
        f_radius: float = FINGERPRINT_RADIUS,
        train_models: bool = True,
        extend_background_samples: bool = False,
        cofactor_augmentation: bool = False,
        threads: int = 1,
    """
    # iterate over each subset
    for subset in os.listdir(directory):
        print("Build fingerprint table for", subset)
        fDB = FingerprintDB()

        subset_path = directory / Path(str(subset))

        train_dir = subset_path / Path("train")
        db_path = subset_path / Path("db")

        fingerprint_db_path = db_path / Path("fingerprintDB")
        os.makedirs(fingerprint_db_path, exist_ok=True)

        out_file = fingerprint_db_path / Path(name + "_fpr" + str(f_radius) + ".tsv")

        fDB.build_db(
            input_directory=train_dir,
            f_radius=f_radius,
            train_models=False,
            extend_background_samples=True,
            cofactor_augmentation=True,
            threads=jobs,
        )
        fDB.save_fingerprint_tsv(out_file)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("dir", type=Path)
    parser.add_argument("--f_radius", type=float, default=5.0)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--name", type=str, default="db")
    args = parser.parse_args()

    build_fingerprint_base(args.dir, args.f_radius, args.jobs, args.name)

"""
creates the fingerprint.tsv in the correct folder
"""

import os
from pathlib import Path
import argparse
from apo2holo.database.fingerprint.build import build_and_safe_fingerprint_tsv


def build_fingerprint_base(
    directory: Path,
    f_radius: float,
    jobs: int,
    name: str,
    fingerprint_types: list,
    extend_background_samples: bool,
    cofactor_augmentation: bool,
    min_dist_art_samples : float,
    n_art_samples : int,
    strategy : str,
    aug_radius : float,
    n_augs : int
) -> str:
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
  
        subset_path = directory / Path(str(subset))

        train_dir = subset_path / Path("train")
        db_path = subset_path / Path("db")

        fingerprint_db_path = db_path / Path("fingerprintDB")
        os.makedirs(fingerprint_db_path, exist_ok=True)

        out_file = fingerprint_db_path / Path(name + "_fpr" + str(f_radius) + ".tsv")

        build_and_safe_fingerprint_tsv(
            input_directory=train_dir,
            out_file=out_file,
            f_radius=f_radius,
            extend_background_samples=extend_background_samples,
            cofactor_augmentation=cofactor_augmentation,
            threads=jobs,
            fingerprint_types=fingerprint_types,
            min_dist_art_samples=min_dist_art_samples,
            n_art_samples=n_art_samples,
            strategy=strategy,
            aug_radius=aug_radius,
            n_augs=n_augs
        )

    return name + "_fpr" + str(f_radius) + ".tsv"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("dir", type=Path)
    parser.add_argument("--f_radius", type=float, default=5.0)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--name", type=str, default="db")
    args = parser.parse_args()

    build_fingerprint_base(args.dir, args.f_radius, args.jobs, args.name)

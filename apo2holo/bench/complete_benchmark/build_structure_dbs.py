"""
create structure databases
"""
import os
from pathlib import Path
import argparse
from apo2holo.build import build_structure_db

def build_structure_dbs(directory : Path):
    """
    """
    # iterate over each subset
    for subset in os.listdir(directory):
        print("Build structure database for", subset)
        subset_path = directory / Path(str(subset))

        train_dir = subset_path / Path('train')

        db_path = subset_path / Path('db')
        os.makedirs(db_path, exist_ok=True)

        chain_dir = db_path / Path('single_chains')
        os.makedirs(chain_dir, exist_ok=True)
        structure_db = db_path / Path('structureDB')
        os.makedirs(structure_db, exist_ok=True)
        structure_db = structure_db / Path('structureDB')

        build_structure_db(train_dir, chain_dir, structure_db, jobs=8)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("dir", type=Path)

    args = parser.parse_args()

    build_structure_dbs(args.dir)
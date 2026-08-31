"""
Builds the databases in order to run the script
"""

import argparse
import os
from src.io.writers.printl import printl
from pathlib import Path
from src.database.structure.StructureDB import StructureDB
from src.database.FingerprintDB import FingerprintDB


def parse_args(args):
    """ """
    input_structure_dir = args.input_structure_dir
    if not os.path.exists(input_structure_dir):
        raise ValueError(f"{input_structure_dir} is no valid path")

    output_structure_dir = args.output_structure_dir
    if not os.path.exists(output_structure_dir):
        printl("creating output folder for chains")
        os.makedirs(output_structure_dir, exist_ok=True)

    structure_db_path = args.structure_db_path
    if structure_db_path:
        if not os.path.exists(structure_db_path):
            printl(f"Creating parent folder {structure_db_path}")
            os.makedirs(structure_db_path, exist_ok=True)
        structure_db_path = os.path.join(structure_db_path, "structureDB")

    return (
        input_structure_dir,
        output_structure_dir,
        structure_db_path,
        args.jobs,
    )


def build(
    input_structure_dir: Path,
    output_structure_dir: Path,
    structure_db_path: Path,
    jobs: int,
):
    """ """
    # build structure db
    build_structure_db(
        input_structure_dir, output_structure_dir, structure_db_path, jobs
    )


def build_structure_db(
    raw_input_dir: Path, chain_dir: Path, structure_db_path: Path, jobs: int
):
    """ """
    structureDB = StructureDB()

    structureDB.build_db(raw_input_dir, chain_dir, structure_db_path, threads=jobs)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "input_structure_dir", help="Path to the folder containing the raw structures."
    )
    parser.add_argument(
        "output_structure_dir", help="Path to the folder for the extracted chains"
    )
    parser.add_argument(
        "--structure-db-path",
        type=str,
        help="Folder for the Structure DB",
        default=None,
    )
    parser.add_argument("--jobs", type=int, help="Amount of cores", default=1)
    args = parser.parse_args()

    (
        input_structure_dir,
        output_structure_dir,
        structure_db_path,
        fingerprint_db_path,
        jobs,
    ) = parse_args(args)
    build(
        input_structure_dir,
        output_structure_dir,
        structure_db_path,
        jobs,
    )

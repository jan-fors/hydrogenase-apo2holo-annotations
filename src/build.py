"""
Builds the databases in order to run the script
"""
import argparse
import os
from src.io.printl import printl
from pathlib import Path
from src.database.StructureDB import StructureDB
from src.database.FingerprintDB import FingerprintDB

def parse_args(args):
    """
    """
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

    fingerprint_db_path = args.fingerprint_db_path
    if fingerprint_db_path:
        if not os.path.exists(fingerprint_db_path):
            printl(f"Creating parent folder {structure_db_path}")
            os.makedirs(fingerprint_db_path, exist_ok=True)

    return input_structure_dir, output_structure_dir, structure_db_path, fingerprint_db_path

def build(input_structure_dir : Path, output_structure_dir : Path, structure_db_path : Path, fingerprint_db_path : Path):
    """
    """
    # build structure db
    #build_structure_db(input_structure_dir, output_structure_dir, structure_db_path)

    # build fingerprint db
    build_fingerprint_db(input_structure_dir, fingerprint_db_path)

def build_structure_db(raw_input_dir : Path, chain_dir : Path, structure_db_path : Path):
    """
    """
    structureDB = StructureDB()

    structureDB.build_db(raw_input_dir, chain_dir, structure_db_path)
    

def build_fingerprint_db(raw_input_dir : Path, fingerprint_db_path : Path):
    """
    """
    fingerprintDB = FingerprintDB(fingerprint_db_path)

    fingerprintDB.build_db(raw_input_dir)

    fingerprintDB.save()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_structure_dir", help="Path to the folder containing the raw structures.")
    parser.add_argument("output_structure_dir", help="Path to the folder for the extracted chains")
    parser.add_argument("--structure-db-path", type=str, help="Folder for the Structure DB", default=None)
    parser.add_argument("--fingerprint-db-path", type=str, help="Path to the fingerprint DB", default=None)
    args = parser.parse_args()

    input_structure_dir, output_structure_dir, structure_db_path, fingerprint_db_path = parse_args(args)
    build(input_structure_dir, output_structure_dir, structure_db_path, fingerprint_db_path)
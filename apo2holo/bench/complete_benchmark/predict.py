"""
Given a chosen model each test data sample is predicted and results are saved
"""

import argparse
import os
from pathlib import Path
import subprocess
from apo2holo.main import main
from datetime import datetime
from typing import Literal


def predict(
    runname : str,
    directory: Path,
    as_model: Path,
    fes_type_model: Path,
    fes_pocket_model: Path,
    as_search_type: Literal["mc", "sum"],
    fes_pocket_search_type : Literal["mc", "sum"],
    fes_type_search_type : Literal["mc", "sum"],
    nn_clustering_radius : float,
    fident_threshold : float,
    bits_threshold : float,
):
    
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

    for subset in os.listdir(directory):
        print("annotate", subset, "...")
        subset_path = directory / Path(str(subset))
        db_dir = subset_path / Path("db")
        structure_db_path = db_dir / Path("structureDB") / Path("structureDB")
        fingerprint_db_dir_path = db_dir / Path("fingerprintDB")
        chain_dir = db_dir/Path("single_chains")
        test_dir = subset_path / Path("test")

        as_model_absolut_path = fingerprint_db_dir_path / as_model
        
        fes_type_model_absolut_path = fingerprint_db_dir_path / fes_type_model

        fes_pocket_model_absolut_path = fingerprint_db_dir_path / fes_pocket_model

        pred_dir = subset_path / Path("out") / Path(timestamp + "_" + runname)

        for test_structure in os.listdir(test_dir):
            try:
                print("Predicting test structure", test_structure)
                test_structure_path = test_dir / Path(test_structure)
                output = pred_dir / Path(test_structure.split(".")[0])
                os.makedirs(output, exist_ok=True)
                main(
                    input_structure_path=test_structure_path,
                    output_dir=output,
                    structure_db_path=structure_db_path,
                    chain_dir_path=chain_dir,
                    nn_clustering_radius = nn_clustering_radius,
                    fident_threshold = fident_threshold,
                    bits_threshold = bits_threshold,
                    as_model_path = as_model_absolut_path,
                    as_search_type = as_search_type,
                    fes_pocket_model_path = fes_pocket_model_absolut_path,
                    fes_pocket_search_type = fes_pocket_search_type,
                    fes_type_model_path = fes_type_model_absolut_path,
                    fes_type_search_type = fes_type_search_type,
                    plot = True,
                    boltz = False
                )
            except (FileNotFoundError, subprocess.CalledProcessError) as e:
                print(e)

    return Path("out") / Path(timestamp + "_" + runname)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("dir", type=Path)
    parser.add_argument("--f_radius", type=float, default=5.0)
    parser.add_argument(
        "--as_model",
        type=Path,
        default=None,
        help="Path inside the fingeprintDB directory to the model that should be used",
    )
    parser.add_argument(
        "--fes_type_model",
        type=Path,
        default=None,
        help="Path inside the fingeprintDB directory to the model that should be used",
    )
    parser.add_argument(
        "--fes_pocket_model",
        type=Path,
        default=None,
        help="Path inside the fingeprintDB directory to the model that should be used",
    )
    parser.add_argument(
        "--model_dir",
        type=Path,
        default=None,
        help="Directory containing multiple model files",
    )
    parser.add_argument("--search_type", type=str, choices=["sum", "mc"], default="sum")
    parser.add_argument("--tmp", type=Path, default=Path("tmp"))

    args = parser.parse_args()

    predict(
        args.dir, args.f_radius, args.as_model, args.fes_type_model, args.fes_pocket_model, args.model_dir, args.search_type, args.tmp
    )

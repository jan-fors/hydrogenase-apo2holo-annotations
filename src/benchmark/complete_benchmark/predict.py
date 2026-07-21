"""
Given a chosen model each test data sample is predicted and results are saved
"""

import argparse
import os
from pathlib import Path
import subprocess
from src.main import main
from datetime import datetime


def predict(
    directory: Path,
    f_radius: float,
    as_model: Path,
    fes_type_model: Path,
    fes_pocket_model: Path,
    model_dir: Path,
    search_type: str,
    tmp: Path,
):
    """
    if a model dir is provided the folder structure should be kept:
        !!! MODEL DIR HAS TO EXIST LIKE THIS IN EVERY SUBSET !!!
        -> input: lrc/R5.0/**models, then output: out/lrc/R5.0/*model*/**results
        - *model* is then: model__<MODEL NUMBER>__f<F-Radius>__<SEARCHTYPE>

    if a single model is provided the folder structure should be:
        - out/<TIMESTAMP>_<MODEL_TYPE>_<SEARCH_TYPE>_<F_RADIUS>_<MODEL_NAME>/**results

    """
    # check whether model dir or model has to be predicted
    if model_dir != None:  # predict dir
        # if not os.path.exists(model_dir):
        #     raise ValueError("Model dir path does not exist")
        # create output dir path
        out_path_per_subset = Path("out") / model_dir
        # iterate over subsets
        for subset in os.listdir(directory):
            print(f"[{datetime.now()}]", "Annotate", subset, "...")

            subset_path = directory / Path(str(subset))
            db_dir = subset_path / Path("db")
            structure_db_path = db_dir / Path("structureDB") / Path("structureDB")
            fingerprint_db_dir_path = db_dir / Path("fingerprintDB")

            test_dir = subset_path / Path("test")

            output_dir = subset_path / out_path_per_subset
            if not os.path.exists(output_dir):
                os.makedirs(output_dir, exist_ok=True)

            model_dir_path = fingerprint_db_dir_path / model_dir
            if not os.path.exists(model_dir_path):
                raise ValueError("model dir does not exists")

            # iterate over models
            for model in os.listdir(model_dir_path):
                model_name = model.split(".")[0]

                run_name = model_name + "__" + str(f_radius) + "__" + search_type

                out_path = output_dir / Path(run_name)
                os.makedirs(out_path, exist_ok=True)

                for test_structure in os.listdir(test_dir):
                    try:
                        print("-> Annotating Structure", test_structure)
                        test_structure_path = test_dir / Path(test_structure)
                        output = out_path / Path(test_structure.split(".")[0])
                        os.makedirs(output, exist_ok=True)
                        main(
                            input_path=test_structure_path,
                            out=output,
                            tmp=tmp,
                            boltz=False,
                            plot=True,
                            structure_db_path=structure_db_path,
                            fingerprint_db_path=fingerprint_db_dir_path,
                            search_type=search_type,
                            model=model_dir_path / Path(model),
                            f_radius=f_radius,
                        )
                    except (FileNotFoundError, subprocess.CalledProcessError) as e:
                        print(e)

    elif as_model != None and fes_type_model != None and fes_pocket_model != None:  # predict with a single model
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        model_name = fes_type_model.stem
        # runname = timestamp+ "__"+ model.parent.name +"__"+
        runname = str(timestamp) + "__" + model_name + "__" + search_type + "__" + str(f_radius)

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

            pred_dir = subset_path / Path("out") / Path(runname)

            for test_structure in os.listdir(test_dir):
                try:
                    print("Predicting test structure", test_structure)
                    test_structure_path = test_dir / Path(test_structure)
                    output = pred_dir / Path(test_structure.split(".")[0])
                    os.makedirs(output, exist_ok=True)
                    main(
                        input_path=test_structure_path,
                        out=output,
                        tmp=tmp,
                        boltz=False,
                        plot=True,
                        structure_db_path=structure_db_path,
                        fingerprint_db_path=fingerprint_db_dir_path,
                        search_type=search_type,
                        chain_dir=chain_dir,
                        active_site_model=as_model_absolut_path,
                        fes_type_model=fes_type_model_absolut_path,
                        fes_pocket_model=fes_pocket_model_absolut_path,
                        f_radius=f_radius,
                    )
                except (FileNotFoundError, subprocess.CalledProcessError) as e:
                    print(e)


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

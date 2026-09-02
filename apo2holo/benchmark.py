""" """

import argparse
from pathlib import Path
import yaml
import os
from datetime import datetime
from apo2holo.bench.complete_benchmark.create_subsets import create_subsets
from apo2holo.bench.complete_benchmark.build_structure_dbs import build_structure_dbs
from apo2holo.bench.complete_benchmark.build_fingerprintdb_base import build_fingerprint_base
from apo2holo.bench.complete_benchmark.deduplicate_fingeprint_db import (
    deduplicate_fingerprint_db,
)
from apo2holo.bench.complete_benchmark.resample_dbs import resample_dbs
from apo2holo.bench.fingerprint_benchmark.registry import perform_gridsearch
from apo2holo.bench.complete_benchmark.train_models import train_model
from apo2holo.bench.complete_benchmark.predict import predict
from apo2holo.bench.complete_benchmark.read_model_results import read_model_results


def benchmark(bench_config: dict, jobs: int, result_dir : Path):
    """ """
    run_name = bench_config["run_name"]

    if bench_config["use_existing_dataset"]:  # dataset already exists
        # check if path to dataset exists
        if not os.path.exists(bench_config["existing_dataset"]):
            raise ValueError(
                f"Path to existing dataset does not exist: {bench_config['existing_dataset']}"
            )
        bench_directory = bench_config["existing_dataset"]

    else:  # need to make all from scratch
        # 1. create subsets
        now = datetime.now()
        f = now.strftime("%Y_%m_%d_%H_%M")
        bench_directory = bench_config["directory"] / Path(f + "_" + run_name)
        os.makedirs(bench_directory, exist_ok=True)

        k = bench_config["k"]
        min_seq_id = bench_config["min_seq_id"]

        test_data_path = bench_config["test_data"]
        if not os.path.exists(test_data_path):
            raise ValueError(f"Test data path {test_data_path} does not exist.")

        train_data_path = bench_config["training_data"]
        if not os.path.exists(train_data_path):
            raise ValueError(f"Train data path {train_data_path} does not exist.")

        create_subsets(
            test_data_dir=test_data_path,
            training_data_dir=train_data_path,
            k=k,
            min_seq_id=min_seq_id,
            out=bench_directory,
        )

        # 2. build structure dbs
        build_structure_dbs(bench_directory)

    # 3. build fingerprintdbs
    fingerprint_db_config = bench_config["fingerprint_db"]
    fingerprint_types = fingerprint_db_config["types"]
    fingerprint_radius = fingerprint_db_config["radius"]

    do_extend_background_samples = fingerprint_db_config["extend_background_samples"][
        "do"
    ]

    min_dist_art_samples = fingerprint_db_config["extend_background_samples"][
        "min_dist_art_samples"
    ]
    n_art_samples = fingerprint_db_config["extend_background_samples"]["n_art_samples"]
    strategy = fingerprint_db_config["extend_background_samples"]["strategy"]

    do_cofactor_augmentation = fingerprint_db_config["cofactor_augmentation"]["do"]

    aug_radius = fingerprint_db_config["cofactor_augmentation"]["aug_radius"]
    n_augs = fingerprint_db_config["cofactor_augmentation"]["n_augs"]

    fingerprint_db_name = build_fingerprint_base(
        directory=bench_directory,
        f_radius=fingerprint_radius,
        jobs=jobs,
        name=run_name,
        extend_background_samples=do_extend_background_samples,
        cofactor_augmentation=do_cofactor_augmentation,
        fingerprint_types=fingerprint_types,
        min_dist_art_samples=min_dist_art_samples,
        n_art_samples=n_art_samples,
        strategy=strategy,
        aug_radius=aug_radius,
        n_augs=n_augs
    )
   
    # deduplicate
    if fingerprint_db_config["deduplicate"]:
        fingerprint_db_name = deduplicate_fingerprint_db(
            bench_directory, fingerprint_db_name
        )

    # pca

    # resample
    if fingerprint_db_config["resample"]["do"]:
        resample_params = fingerprint_db_config["resample"]["params"]
        undersampling = resample_params["undersampling"]
        top = resample_params["top"]
        fingerprint_db_name = resample_dbs(
            bench_directory, fingerprint_db_name, undersampling, top
        )

    # 4. perform gridsearch on each dataset and choose best params
    model_params = bench_config["models"]
    specific_fingerprint_db_path = (
        bench_directory
        / Path("0")
        / Path("db")
        / Path("fingerprintDB")
        / Path(fingerprint_db_name)
    )

    # active_site
    active_site_params = model_params["active_site"]

    active_site_gridsearch_results = perform_gridsearch(
        model_type=active_site_params["type"],
        data=specific_fingerprint_db_path,
        type="as",
        n_iter=active_site_params["n_iter"],
        jobs=jobs,
    )

    # fes pocket
    fes_pocket_params = model_params["fes"]["pocket"]

    fes_pocket_gridsearch_results = perform_gridsearch(
        model_type=fes_pocket_params["type"],
        data=specific_fingerprint_db_path,
        type="fes_pocket",
        n_iter=fes_pocket_params["n_iter"],
        jobs=jobs,
    )

    # fes type
    fes_type_params = model_params["fes"]["type"]

    fes_type_gridsearch_results = perform_gridsearch(
        model_type=fes_type_params["type"],
        data=specific_fingerprint_db_path,
        type="fes",
        n_iter=fes_type_params["n_iter"],
        jobs=jobs,
    )

    # 5. train models with params
    model_dir_name = "model"
    # train as model
    as_model = train_model(
        directory=bench_directory,
        finerprinttsv_name=fingerprint_db_name,
        model_name=run_name,
        model_type=active_site_params["type"],
        model_dir_name=model_dir_name,
        training_params=active_site_gridsearch_results["best_params"],
        pred_type="as",
        fingerprint_type=fingerprint_types,
        f_radius=fingerprint_radius,
    )

    # train fes pocket model
    fes_pocket_model = train_model(
        directory=bench_directory,
        finerprinttsv_name=fingerprint_db_name,
        model_name=run_name,
        model_dir_name=model_dir_name,
        model_type=fes_pocket_params["type"],
        training_params=fes_pocket_gridsearch_results["best_params"],
        pred_type="fes_pocket",
        fingerprint_type=fingerprint_types,
        f_radius=fingerprint_radius,
    )

    # train fes type model
    fes_type_model = train_model(
        directory=bench_directory,
        finerprinttsv_name=fingerprint_db_name,
        model_name=run_name,
        model_dir_name=model_dir_name,
        model_type=fes_type_params["type"],
        training_params=fes_type_gridsearch_results["best_params"],
        pred_type="fes",
        fingerprint_type=fingerprint_types,
        f_radius=fingerprint_radius,
    )

    # 6. predict
    out_dir = predict(
        runname=run_name,
        directory=bench_directory,
        as_model=as_model,
        fes_pocket_model=fes_pocket_model,
        fes_type_model=fes_type_model,
        as_search_type=bench_config["predict"]["search_type"],
        fes_pocket_search_type=bench_config["predict"]["search_type"],
        fes_type_search_type=bench_config["predict"]["search_type"],
        nn_clustering_radius=bench_config["structure_db"]["nn_clustering_radius"],
        fident_threshold=bench_config["structure_db"]["FIDENT_THRESHOLD"],
        bits_threshold=bench_config["structure_db"]["BITS_THRESHOLD"],
    )

    # 7. analyze results
    # read model results
    read_model_results(bench_directory, out_dir, result_dir)

    # score model results

    # save best model for inference if SAVE


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, help="Path to the benchmark config")
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--result_dir", type=Path, default="benchmarking")
    args = parser.parse_args()

    if not os.path.exists(args.config):
        raise ValueError(f"{args.config} does not exist.")

    with open(args.config, "r") as file:
        config = yaml.safe_load(file)

    if not os.path.exists(args.result_dir):
        os.makedirs(args.result_dir, exist_ok=True)

    benchmark(config, args.jobs, args.result_dir)

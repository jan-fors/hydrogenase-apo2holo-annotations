import argparse
from pathlib import Path
import os
from src.io.readers.read_config import read_config

from src.updated_main import main


def cli(args):
    """ """
    input_file = args.input
    if not os.path.exists(input_file):
        raise ValueError(f"Input file: {input_file} does not exist")

    config = args.config
    if not os.path.exists(config):
        raise ValueError(f"Config file: {config} does not exist")

    # create specific output folder
    output_dir = args.output_dir
    specific_output_dir = output_dir / Path(input_file.stem)
    if not os.path.exists(specific_output_dir):
        os.makedirs(specific_output_dir, exist_ok=True)

    config_data = read_config(config)

    structure_db_path = config_data["structure_db"]["path"]
    chain_dir_path = config_data["structure_db"]["chain_dir"]
    nn_clustering_radius = config_data["structure_db"]["nn_clustering_radius"]
    fident_threshold = config_data["structure_db"]["FIDENT_THRESHOLD"]
    bits_threshold = config_data["structure_db"]["BITS_THRESHOLD"]

    as_model_path = config_data["models"]["active_site"]["path"]
    as_search_type = config_data["models"]["active_site"]["search_type"]

    fes_pocket_model_path = config_data["models"]["fes"]["pocket"]["path"]
    fes_pocket_search_type = config_data["models"]["fes"]["pocket"]["search_type"]

    fes_type_model_path = config_data["models"]["fes"]["type"]["path"]
    fes_type_search_type = config_data["models"]["fes"]["type"]["search_type"]


    main(
        input_structure_path=input_file,
        output_dir=specific_output_dir,
        structure_db_path=structure_db_path,
        chain_dir_path=chain_dir_path,
        nn_clustering_radius=nn_clustering_radius,
        fident_threshold=fident_threshold,
        bits_threshold=bits_threshold,
        as_model_path=as_model_path,
        as_search_type=as_search_type,
        fes_pocket_model_path=fes_pocket_model_path,
        fes_pocket_search_type=fes_pocket_search_type,
        fes_type_model_path=fes_type_model_path,
        fes_type_search_type=fes_type_search_type,
        plot=args.plot,
        boltz=args.boltz
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("input", type=Path, help="Input structure")
    parser.add_argument(
        "-c",
        "--config",
        type=Path,
        help="Path to the config file",
        default="config/inference/config.yml",
    )
    parser.add_argument(
        "-o",
        "--output_dir",
        type=Path,
        help="Path to the output directory",
        default=".",
    )
    parser.add_argument("--plot", action="store_true")
    parser.add_argument("-boltz", action="store_true")

    args = parser.parse_args()

    cli(args)

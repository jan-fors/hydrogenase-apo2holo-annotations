import argparse
from pathlib import Path
import os
from src.main import main
from src.utils.constants import (
    STRUCTURE_DB,
    FINGERPRINT_DB,
    SEARCH_TYPE,
    STRUCTURE_DIR,
    FINGERPRINT_RADIUS,
    FES_POCKET_MODEL,
    FES_TYPE_MODEL,
    ACTIVE_SITE_MODEL,
)


def _verify_inputs(args):
    """"""
    # TODO check if db exists

    # TODO check if input file exists and has valid format

    pass


def _extract_args(args):
    """"""
    if args.active_site_model != None:
        active_site_model = Path(args.active_site_model)
    else:
        active_site_model = None

    if args.fes_type_model != None:
        fes_type_model = Path(args.fes_type_model)
    else:
        fes_type_model = None

    if args.fes_pocket_model != None:
        fes_pocket_model = Path(args.fes_pocket_model)
    else:
        fes_pocket_model = None

    return (
        Path(args.input_path),
        args.output,
        args.output_dir,
        args.tmp,
        args.boltz,
        args.plot,
        Path(args.structure_db_path),
        Path(args.fingerprint_db_path),
        args.search_type,
        active_site_model,
        fes_type_model,
        fes_pocket_model,
        args.f_radius,
        args.chain_dir,
    )


def cli(args):
    """ """

    _verify_inputs(args)

    # extract inputs
    (
        input_path,
        output,
        output_dir,
        tmp,
        boltz,
        plot,
        structure_db_path,
        fingerprint_db_path,
        search_type,
        active_site_model,
        fes_type_model,
        fes_pocket_model,
        f_radius,
        chain_dir,
    ) = _extract_args(args)

    # create folders if necessary
    out = Path(os.path.join(output_dir, output))
    os.makedirs(out, exist_ok=True)

    main(
        input_path=input_path,
        out=out,
        tmp=tmp,
        boltz=boltz,
        plot=plot,
        structure_db_path=structure_db_path,
        fingerprint_db_path=fingerprint_db_path,
        search_type=search_type,
        active_site_model=active_site_model,
        fes_type_model=fes_type_model,
        fes_pocket_model=fes_pocket_model,
        f_radius=f_radius,
        chain_dir=chain_dir,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path", help="Path to the input structure.")
    parser.add_argument("output", help="Name of the output.")
    parser.add_argument(
        "--output_dir",
        "-o",
        default=".",
        help="Specify the output directory, default .",
    )
    parser.add_argument("--tmp", default="tmp", help="Path to the tmp folder.")
    parser.add_argument(
        "--boltz", action="store_true", help="Create Boltz input yaml file."
    )
    parser.add_argument(
        "--plot", action="store_true", help="Create a plot of the results."
    )
    parser.add_argument(
        "--result-table", default=None, help="Path to a result table for batch runs."
    )
    parser.add_argument(
        "--structure-db-path",
        default=STRUCTURE_DB,
        help="Path to the structure database",
    )
    parser.add_argument(
        "--fingerprint-db-path",
        default=FINGERPRINT_DB,
        help="Path to the fingerprint database",
    )
    parser.add_argument("--search_type", type=str, default=SEARCH_TYPE)
    parser.add_argument("--active_site_model", type=Path, default=ACTIVE_SITE_MODEL)
    parser.add_argument("--fes_type_model", type=Path, default=FES_TYPE_MODEL)
    parser.add_argument("--fes_pocket_model", type=Path, default=FES_POCKET_MODEL)
    parser.add_argument("--f_radius", type=float, default=FINGERPRINT_RADIUS)
    parser.add_argument("--chain_dir", type=Path, default=STRUCTURE_DIR)

    args = parser.parse_args()

    cli(args)

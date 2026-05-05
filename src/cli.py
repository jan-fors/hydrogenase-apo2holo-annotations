import argparse
from pathlib import Path
import os
from src.main import main
from src.utils.constants import (
    STRUCTURE_DB,
    FINGERPRINT_DB,
    SEARCH_TYPE
)

def _verify_inputs(args):
    """"""
    #TODO check if db exists

    #TODO check if input file exists and has valid format

    pass

def _extract_args(args):
    """"""
    return Path(args.input_path), args.output, args.output_dir, args.tmp, args.boltz, args.plot, args.result_table, Path(args.structure_db_path), Path(args.fingerprint_db_path), args.search_type

def cli(args):
    """
    """

    #    
    _verify_inputs(args)

    # extract inputs
    input_path, output, output_dir, tmp, boltz, plot, result_table, structure_db_path, fingerprint_db_path, search_type = _extract_args(args)

    # create folders if necessary
    out = Path(os.path.join(output_dir, output))
    os.makedirs(out, exist_ok=True)

    main(input_path=input_path, 
         out=out, 
         tmp=tmp, 
         boltz=boltz, 
         plot=plot,
         result_table_path=result_table,
         structure_db_path=structure_db_path,
         fingerprint_db_path=fingerprint_db_path,
         search_type=search_type)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path", help="Path to the input structure.")
    parser.add_argument("output", help="Name of the output.")
    parser.add_argument("--output_dir", "-o", default=".", help="Specify the output directory, default .")
    parser.add_argument("--tmp", default="tmp", help="Path to the tmp folder.")
    parser.add_argument("--boltz", action="store_true", help="Create Boltz input yaml file.")
    parser.add_argument("--plot", action="store_true", help="Create a plot of the results.")
    parser.add_argument("--result-table", default=None, help="Path to a result table for batch runs.")
    parser.add_argument("--structure-db-path", default=STRUCTURE_DB, help="Path to the structure database")
    parser.add_argument("--fingerprint-db-path", default=FINGERPRINT_DB, help="Path to the fingerprint database")
    parser.add_argument("--search-type", default=SEARCH_TYPE)
    # add db params and model param

    args = parser.parse_args()

    cli(args)
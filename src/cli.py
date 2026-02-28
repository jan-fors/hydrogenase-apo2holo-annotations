import argparse
from src.main import main

def _verify_inputs(args):
    """"""
    #TODO check if db exists

    #TODO check if input file exists and has valid format

    pass

def _extract_args(args):
    """"""
    return args.input_path, args.output, args.output_dir, args.tmp, args.boltz

def cli(args):
    """
    """
    
    _verify_inputs(args)

    input_path, output, output_dir, tmp, boltz = _extract_args(args)

    main(input_path, output, output_dir, tmp, boltz)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input_path", help="Path to the input structure.")
    parser.add_argument("output", help="Name of the output.")
    parser.add_argument("--output_dir", "-o", default=".", help="Specify the output directory, default .")
    parser.add_argument("--tmp", default="tmp", help="Path to the tmp folder.")
    parser.add_argument("--boltz", action="store_true", help="Create Boltz input yaml file.")
    args = parser.parse_args()

    cli(args)
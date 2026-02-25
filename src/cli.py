import argparse
from src.main import main

def _verify_inputs(args):
    """"""
    #TODO check if db exists

    #TODO check if input file exists and has valid format

    pass

def cli(args):
    """
    """
    
    _verify_inputs(args)

    
    main(args)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("input", help="Path to the input structure.")
    parser.add_argument("output", help="Name of the output.")
    parser.add_argument("--output-dir", "-o", default=".", help="Specify the output directory, default .")
    parser.add_argument("--database", default="/db", help="Path to the database")
    parser.add_argument("--boltz", action="store_true", help="Create Boltz input yaml file.")
    args = parser.parse_args()

    cli(args)
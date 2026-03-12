from pathlib import Path
import subprocess
from src.utils.constants import (
    STRUCTURE_DB,
    FOLDSEEK_OUT_FORMAT
)
import os

def search_against_structuredb(input_path : Path, output_path : Path, tmp_path : Path, structure_db_path : str = STRUCTURE_DB):
    """
    foldseek easy-multimersearch example/1tim.pdb.gz example/8tim.pdb.gz result tmpFolder
    """
    result_name = os.path.basename(input_path).split(".")[0] + "_ms_res"
    output_path = os.path.join(str(output_path), result_name)

    cmd = ["foldseek", "easy-search", str(input_path), structure_db_path, str(output_path), str(tmp_path), "--format-output", FOLDSEEK_OUT_FORMAT] #TODO add other arguments

    subprocess.run(cmd, stderr=subprocess.PIPE, text=True, check=True)

    return output_path
import os
from pathlib import Path
import subprocess

def extract_chain(input_pdb: str, output_dir: str, chain: str):
    """ """
    structure_name = os.path.basename(input_pdb).split(".")[0] + f"_{chain}.pdb"
    outpath = Path(os.path.join(output_dir, structure_name))
    cmd = ["pdb_selchain", f"-{chain}", str(input_pdb)]

    with outpath.open("w") as f:
        subprocess.run(cmd, stdout=f, stderr=subprocess.PIPE, text=True, check=True)

    return outpath
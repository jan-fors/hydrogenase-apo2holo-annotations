from src.utils.constants import (
    STRUCTURE_DIR
)
from pathlib import Path
import os


def get_structure_path(name : str) -> Path:
    """"""
    structure_path = os.path.join(STRUCTURE_DIR, name +".pdb")
    if not os.path.exists(structure_path):
        print(f"{structure_path} does not exist.")
        return None
    return Path(structure_path)
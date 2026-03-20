from src.utils.constants import (
    COFACTOR_BLACKLIST
)
from pathlib import Path
from src.io.printl import printl
import os

def apply_blacklist_build(cofactors: dict[str, list[tuple[str, float, float, float]]]) -> dict[str, list[tuple[str, float, float, float]]]:
    res = {}
    black_listed = {}
    for key in cofactors.keys():
        # TODO change later but for testing initially following is okay
        if cofactors[key]["res_name"] in COFACTOR_BLACKLIST:
            black_listed[key] = cofactors[key]
            continue
        res[key] = cofactors[key]
    return res, black_listed

def apply_blacklist(cofactors: dict[str, list[tuple[str, float, float, float]]]) -> dict[str, list[tuple[str, float, float, float]]]:
    res = {}
    
    for key in cofactors.keys():
        # TODO change later but for testing initially following is okay
        if key in COFACTOR_BLACKLIST:
            continue
        res[key] = cofactors[key]
    return res

def apply_blacklist_to_input_structures(structure_path: str):
    """ """
    pdb_path = Path(structure_path)
    tmp_path = pdb_path.with_suffix(structure_path.suffix + ".tmp")

    removed = 0

    with pdb_path.open("r") as fin, tmp_path.open("w") as fout:
        for line in fin:
            record = line[0:6].strip()
            if record == "HETATM":
                resname = line[17:20].strip().upper()
                if resname in COFACTOR_BLACKLIST:
                    removed += 1
                    continue
            fout.write(line)

    # Atomar ersetzen (sehr wichtig!)
    os.replace(tmp_path, pdb_path)

    printl(f"[remove_blacklisted_hetatm_inplace] removed_hetatm_lines={removed}")
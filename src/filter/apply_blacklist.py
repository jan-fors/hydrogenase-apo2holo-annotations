from src.utils.constants import (
    COFACTOR_BLACKLIST
)

def apply_blacklist_build(cofactors: dict[str, list[tuple[str, float, float, float]]]) -> dict[str, list[tuple[str, float, float, float]]]:
    res = {}
    
    for key in cofactors.keys():
        # TODO change later but for testing initially following is okay
        if cofactors[key]["res_name"] in COFACTOR_BLACKLIST:
            continue
        res[key] = cofactors[key]
    return res

def apply_blacklist(cofactors: dict[str, list[tuple[str, float, float, float]]]) -> dict[str, list[tuple[str, float, float, float]]]:
    res = {}
    
    for key in cofactors.keys():
        # TODO change later but for testing initially following is okay
        if key in COFACTOR_BLACKLIST:
            continue
        res[key] = cofactors[key]
    return res
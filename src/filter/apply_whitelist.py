from src.utils.constants import (
    COFACTOR_WHITELIST
)

def apply_whitelist(cofactors: dict[str, list[tuple[str, float, float, float]]]) -> dict[str, list[tuple[str, float, float, float]]]:
    res = {}
    
    for key in cofactors.keys():
        if key in COFACTOR_WHITELIST:
            res[key] = cofactors[key]
    return res

def apply_whitelist_build(cofactors: dict[str, list[tuple[str, float, float, float]]]) -> dict[str, list[tuple[str, float, float, float]]]:
    res = {}
    not_whitelisted = {}

    for key in cofactors.keys():
        if cofactors[key]["res_name"] in COFACTOR_WHITELIST:
            res[key] = cofactors[key]
        else:
            not_whitelisted[key] = cofactors[key]
    return res, not_whitelisted
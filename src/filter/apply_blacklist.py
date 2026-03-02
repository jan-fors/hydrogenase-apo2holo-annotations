from src.utils.constants import (
    COFACTOR_BLACKLIST
)

def apply_blacklist(cofactors: dict[str, list[tuple[str, float, float, float]]]) -> dict[str, list[tuple[str, float, float, float]]]:
    res = {}

    for key in cofactors.keys():
        if key in COFACTOR_BLACKLIST:
            continue
        res[key] = cofactors[key]
    return res
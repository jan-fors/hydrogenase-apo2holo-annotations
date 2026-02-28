from pathlib import Path
from typing import List

def extract_cofactors(pdb_path: str) -> dict[str, list[tuple[str, float, float, float]]]:
    """
    Returns:
        {
            "RESNAME": [(element, x, y, z), ...],
            ...
        }
    """

    het_groups = {}

    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("HETATM"):

                resname = line[17:20].strip()
                element = line[76:78].strip()

                # Fallback: falls Element-Spalte leer ist
                if not element:
                    element = line[12:16].strip()[0]

                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])

                if resname not in het_groups:
                    het_groups[resname] = []

                het_groups[resname].append((element, x, y, z))

    return het_groups
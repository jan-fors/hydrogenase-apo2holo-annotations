from typing import List

def get_chains(input_pdb: str) -> List[str]:
    """ """
    chains = set()
    with open(input_pdb, "r") as f:
        for line in f:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 22:
                c = line[21].strip()
                if c:
                    chains.add(c)
    return sorted(chains)
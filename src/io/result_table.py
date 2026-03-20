import csv

def write_to_result_table(data : dict, path : str):
    """
    Write specified data into the result table
    """
    new_line = {
        "active_site": data["active_site"]["smiles"],
        "proximal": data["proximal"]["smiles"],
        "medial": data["medial"]["smiles"],
        "distal": data["distal"]["smiles"]
    }

    with open(path, "a", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["active_site", "proximal", "medial", "distal"])
        writer.writerow(new_line)
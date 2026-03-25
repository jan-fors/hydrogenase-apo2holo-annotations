import csv

def write_to_result_table(data : dict, path : str):
    """
    Write specified data into the result table
    """
    if "proximal" in data.keys():
        new_line = {
            "active_site": data["active_site"]["formula"],
            "proximal": data["proximal"]["formula"],
            "medial": data["medial"]["formula"],
            "distal": data["distal"]["formula"]
        }

        with open(path, "a", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=["active_site", "proximal", "medial", "distal"])
            writer.writerow(new_line)

    else:
        new_line = {}
        for key in data.keys():
            new_line[key] = data[key]["formula"]

        with open(path, "a", newline="", encoding="utf-8") as f:
            writer = csv.DictWriter(f, fieldnames=data.keys())
            writer.writerow(new_line)


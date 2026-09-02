from pathlib import Path
from typing import List
from apo2holo.utils.formula.smiles import SMILES
import os
import json
from apo2holo.io.writers.printl import printl
from apo2holo.pdb.pdb_handler import identify_cysteines
import numpy as np
from apo2holo.utils.geometric.calculate_geometric_centers import calculate_center
from apo2holo.io.writers.plot_handler import plot_protein

def write_outputs(
    out: Path,
    sorted_list: List,
    structure_path: Path,
    as_pockets: dict,
    as_pocket_key: int,
    as_probabilities: dict,
    fes_pockets: dict,
    prediction_per_pocket: dict,
    best_hits_per_pocket: dict,
    proximal_key: int,
    medial_key: int,
    distal_key: int,
    boltz: bool,
    plot: bool,
):
    """"""
    # raw_results
    output_dict = {}

    # active_site
    output_dict["active_site"] = {
        "prediction": as_probabilities[as_pocket_key],
        "cysteines": identify_cysteines(
            structure_path=structure_path,
            coords=calculate_center(np.array(as_pockets[as_pocket_key])),
            radius=5.0,
        ),
    }

    for k in fes_pockets.keys():
        # print(k)
        if best_hits_per_pocket[k] != "protein":
            cysteines = identify_cysteines(
                structure_path=structure_path,
                coords=calculate_center(np.array(fes_pockets[k])),
                radius=5.0,
            )
            output_dict[k] = {
                "cluster": int(k),
                "predicted_class": best_hits_per_pocket[k],
                "smiles": SMILES[best_hits_per_pocket[k]],
                "prediction": prediction_per_pocket[k],
                "cysteines": cysteines,
            }

    output_list = []
    output_list.append(output_dict["active_site"])

    for i in sorted_list:
        if int(i[0]) in output_dict.keys():
            output_list.append(output_dict[i[0]])

    ## write json output file
    json_path = os.path.join(out, "raw_results.json")
    json_text = json.dumps(output_list, ensure_ascii=False, indent=4)
    printl(f"Writing raw output to {json_path}")
    with open(json_path, "w", encoding="utf-8") as f:
        f.write(json_text)

    # results
    result_dict = {}

    result_dict["active_site"] = (output_dict["active_site"],)
    if proximal_key != None:
        result_dict["proximal"] = (output_dict[proximal_key],)
    if medial_key != None:
        result_dict["medial"] = (output_dict[medial_key],)
    if distal_key != None:
        result_dict["distal"] = output_dict[distal_key]

    if len(result_dict.keys()) > 0:
        json_path = os.path.join(out, "results.json")
        json_text = json.dumps(result_dict, ensure_ascii=False, indent=4)
        printl(f"Writing structured output to {json_path}")
        with open(json_path, "w", encoding="utf-8") as f:
            f.write(json_text)

    # write boltz yaml output file
    if boltz:
        print("YAML CREATION TODO")

    # create output graphic?
    if plot:
        plot_protein(
            structure_path=structure_path,
            pockets=fes_pockets,
            pred_per_pocket=prediction_per_pocket,
            output=out,
        )
import yaml
from pathlib import Path

ID_LIST = ["Z", "Y", "X", "W", "V"]

class FlowList(list):
    pass

def flow_list_representer(dumper, data):
    return dumper.represent_sequence("tag:yaml.org,2002:seq", data, flow_style=True)

yaml.add_representer(FlowList, flow_list_representer)

def _add_sequence(data_dict : dict, chain : str, sequence : str) -> dict:
    """
    """
    data_dict["sequences"].append({"protein":{"id":chain, "sequence":str(sequence), "msa": None}})
    return data_dict

def _add_ligand(data_dict : dict, smiles : str, id : str) -> dict:
    """
    """
    data_dict["sequences"].append({"ligand":{"smiles":smiles, "id":[id]}})
    return data_dict

def _add_constraint(data_dict: dict, binder : str, contacts : list)->dict:
    """
    """
    data_dict["constraints"].append({"pocket": {"binder":binder, "contacts": contacts}})
    return data_dict

def write_boltz_output(sequences : dict, result_dict : dict, out : Path):
    """
    Sequences needs to be a dict with chain_id : sequence.
    result_dict contains the apo2holo results as well as the constraints etc.
    Out is the output dir
    """
    data = {
        "version": 1,
        "sequences":[
        ]
    }
    # add sequences
    for seq in sequences.keys():
        data = _add_sequence(data, seq, sequences[seq])

    # add ligands
    counter = 0
    for cof in result_dict.keys():
        data = _add_ligand(data, result_dict[cof]["smiles"], str(ID_LIST[counter]))
        counter += 1

    counter = 0

    # add constraints
    data["constraints"] = []

    for cof in result_dict.keys():
        cysteines = FlowList([FlowList([c[0], c[1]]) for c in result_dict[cof]["cysteines"]])
        data = _add_constraint(data, ID_LIST[counter], cysteines)
        counter += 1


    with open(out/Path("result.yml"), "w") as f:
        yaml.dump(data, f)


"""
Takes a set and model as input and creates a tsv file containing the results
"""
import argparse
import os
from pathlib import Path
import json
import pandas as pd

def read_model_results(directory : Path, prediction_out_dir : Path, output_dir : Path):
    """
    1. iterate over each subset
    2. get model dir
    3. read all results into dataframe
    4. save dataframe
    """

    res = pd.DataFrame()
    proba_res = pd.DataFrame()

    for subset in os.listdir(directory):
        model_results = directory / Path(subset) / prediction_out_dir

        for hyd in os.listdir(model_results):
            # result classes
            hyd_results = {}
            hyd_results["hyd"] = [hyd]
            results_path = model_results / Path(hyd) / Path("results.json")

            if not os.path.exists(results_path):
                hyd_results["proximal"] = ["n"]
                hyd_results["medial"] = ["n"]
                hyd_results["distal"] = ["n"]
            else:
                with open(results_path, "r") as f:
                    data = json.loads(f.read())

                try:
                    if "proximal" in data.keys():
                        if type(data["proximal"]) == list:
                            hyd_results["proximal"] = [data["proximal"][0]["predicted_class"]]
                        else:
                            hyd_results["proximal"] = [data["proximal"]["predicted_class"]]
                    else:
                        hyd_results["proximal"] = ["n"]
                except:
                    hyd_results["proximal"] = ["n"]

                try:
                    if "medial" in data.keys():
                        if type(data["medial"]) == list:
                            hyd_results["medial"] = [data["medial"][0]["predicted_class"]]
                        else:
                            hyd_results["medial"] = [data["medial"]["predicted_class"]]
                    else:
                        hyd_results["medial"] = ["n"]
                except:
                    hyd_results["medial"] = ["n"]

                try:
                    if "distal" in data.keys():
                        if type(data["distal"]) == list:
                            hyd_results["distal"] = [data["distal"][0]["predicted_class"]]
                        else:
                            hyd_results["distal"] = [data["distal"]["predicted_class"]]
                    else:
                        hyd_results["distal"] = ["n"]
                except:
                    hyd_results["distal"] = ["n"]

            res = pd.concat([res, pd.DataFrame(hyd_results)], ignore_index = True)  

            # result probas
            proba_results = {}
            proba_results["hyd"] = hyd

            if not os.path.exists(results_path):
                proba_results["proximal_3FE4S"] = None
                proba_results["proximal_4FE4S"] = None
                proba_results["proximal_4FE3S"] = None

                proba_results["medial_3FE4S"] = None
                proba_results["medial_4FE4S"] = None
                proba_results["medial_4FE3S"] = None

                proba_results["distal_3FE4S"] = None
                proba_results["distal_4FE4S"] = None
                proba_results["distal_4FE3S"] = None

            else:
                with open(results_path, "r") as f:
                    data = json.loads(f.read())
                try:
                    if "proximal" in data.keys():
                        if type(data["proximal"]) == list:
                            proba_results["proximal_3FE4S"] = [data["proximal"][0]["prediction"]["3FE4S"]]
                            proba_results["proximal_4FE4S"] = [data["proximal"][0]["prediction"]["4FE4S"]]
                            proba_results["proximal_4FE3S"] = [data["proximal"][0]["prediction"]["4FE3S"]]
                        else:
                            proba_results["proximal_3FE4S"] = [data["proximal"]["prediction"]["3FE4S"]]
                            proba_results["proximal_4FE4S"] = [data["proximal"]["prediction"]["4FE4S"]]
                            proba_results["proximal_4FE3S"] = [data["proximal"]["prediction"]["4FE3S"]]
                    else:
                        proba_results["proximal_3FE4S"] = None
                        proba_results["proximal_4FE4S"] = None
                        proba_results["proximal_4FE3S"] = None
                except:
                    proba_results["proximal_3FE4S"] = None
                    proba_results["proximal_4FE4S"] = None
                    proba_results["proximal_4FE3S"] = None

                try:
                    if "medial" in data.keys():
                        if type(data["medial"]) == list:
                            proba_results["medial_3FE4S"] = [data["medial"][0]["prediction"]["3FE4S"]]
                            proba_results["medial_4FE4S"] = [data["medial"][0]["prediction"]["4FE4S"]]
                            proba_results["medial_4FE3S"] = [data["medial"][0]["prediction"]["4FE3S"]]
                        else:    
                            proba_results["medial_3FE4S"] = [data["medial"]["prediction"]["3FE4S"]]
                            proba_results["medial_4FE4S"] = [data["medial"]["prediction"]["4FE4S"]]
                            proba_results["medial_4FE3S"] = [data["medial"]["prediction"]["4FE3S"]]
                    else:
                        proba_results["medial_3FE4S"] = None
                        proba_results["medial_4FE4S"] = None
                        proba_results["medial_4FE3S"] = None
                except:
                    proba_results["medial_3FE4S"] = None
                    proba_results["medial_4FE4S"] = None
                    proba_results["medial_4FE3S"] = None

                try:
                    if "distal" in data.keys():
                        if type(data["distal"]) == list:
                            proba_results["distal_3FE4S"] = [data["distal"][0]["prediction"]["3FE4S"]]
                            proba_results["distal_4FE4S"] = [data["distal"][0]["prediction"]["4FE4S"]]
                            proba_results["distal_4FE3S"] = [data["distal"][0]["prediction"]["4FE3S"]]   
                        else: 
                            proba_results["distal_3FE4S"] = [data["distal"]["prediction"]["3FE4S"]]
                            proba_results["distal_4FE4S"] = [data["distal"]["prediction"]["4FE4S"]]
                            proba_results["distal_4FE3S"] = [data["distal"]["prediction"]["4FE3S"]]
                    else:
                        proba_results["distal_3FE4S"] = None
                        proba_results["distal_4FE4S"] = None
                        proba_results["distal_4FE3S"] = None
                except:
                    proba_results["distal_3FE4S"] = None
                    proba_results["distal_4FE4S"] = None
                    proba_results["distal_4FE3S"] = None

                #print(proba_results)

                proba_res = pd.concat([proba_res, pd.DataFrame(proba_results)], ignore_index = True) 

    out_path =  output_dir / Path(str(prediction_out_dir.stem) + "_hyd_results.tsv")
    res.to_csv(out_path, sep="\t", index=False)

    out_path =  output_dir / Path(str(prediction_out_dir.stem) + "_proba_results.tsv")
    proba_res.to_csv(out_path, sep="\t", index=False)
        


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=Path)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--output_dir", type=Path, default=".")
    args = parser.parse_args()

    read_model_results(args.dir, args.model, args.output_dir)
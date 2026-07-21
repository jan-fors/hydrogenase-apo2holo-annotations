"""
Takes a set and model as input and creates a tsv file containing the results
"""
import argparse
import os
from pathlib import Path
import json
import pandas as pd

def read_model_results(directory : Path, model : Path, output_dir : Path):
    """
    1. iterate over each subset
    2. get model dir
    3. read all results into dataframe
    4. save dataframe
    """

    res = pd.DataFrame()

    for subset in os.listdir(directory):
        model_results = directory / Path(subset) / Path("out") /  model

        if not os.path.exists(model_results):
            print(f"{model_results} does not exist")
            continue
        else:
            for hyd in os.listdir(model_results):
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
                        
                    if "proximal" in data.keys():
                        hyd_results["proximal"] = [data["proximal"][0]["predicted_class"]]
                    else:
                        hyd_results["proximal"] = ["n"]

                    if "medial" in data.keys():
                        hyd_results["medial"] = [data["medial"][0]["predicted_class"]]
                    else:
                        hyd_results["medial"] = ["n"]

                    if "distal" in data.keys():
                        hyd_results["distal"] = [data["distal"]["predicted_class"]]
                    else:
                        hyd_results["distal"] = ["n"]

                res = pd.concat([res, pd.DataFrame(hyd_results)], ignore_index = True)    

    out_path =  output_dir / Path(str(model.stem) + "_results.tsv")
    res.to_csv(out_path, sep="\t", index=False)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=Path)
    parser.add_argument("--model", type=Path, required=True)
    parser.add_argument("--output_dir", type=Path, default=".")
    args = parser.parse_args()

    read_model_results(args.dir, args.model, args.output_dir)
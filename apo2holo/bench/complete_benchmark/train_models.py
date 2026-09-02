"""
Creates a specific model for each training dataset

Create a list of params for the models and generate all of them to check afterwards

Allow to generate multiple models sampling from a distribution of parameters for RandomizedCV
"""
import os
os.environ["OMP_NUM_THREADS"] = "1"
os.environ["OPENBLAS_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"
os.environ["NUMEXPR_NUM_THREADS"] = "1"
os.environ["VECLIB_MAXIMUM_THREADS"] = "1" 

from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import pickle
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from scipy.stats import loguniform, uniform, randint
import datetime
import warnings
from apo2holo.models.architectures.registry import build_model
from apo2holo.io.writers.save_model import save_model
from apo2holo.utils.constants import COFACTOR_BLACKLIST
warnings.filterwarnings("ignore", module="sklearn")



# FUNCTIONS
def printl(text: str):
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{now}] {text}")

def get_pipeline(model):
    """ """
    # define model
    pipe = Pipeline(
        [
            ("scaler", StandardScaler()),
            ("clf", model),
        ]
    )
    return pipe


def prepare_data(data: Path, type : str):
    """ """
    df = pd.read_csv(data, sep="\t")
    printl("Read data...")

    if "id" in df.columns:
        df.drop("id", axis=1, inplace=True)
    if "res_name" in df.columns:
        df.drop("res_name", axis=1, inplace=True)
    if "smiles" in df.columns:
        df.drop("smiles", axis=1, inplace=True)
    if "type" in df.columns:
        df = df.drop(columns="type")
    if "aug_dist" in df.columns:
        df = df.drop(columns="aug_dist")
    if "structure" in df.columns:
        df = df.drop(columns="structure")
    if "fingerprint" in df.columns:
        df.drop("fingerprint", axis=1, inplace=True)
    df = df[~df["formula"].isin(COFACTOR_BLACKLIST)].copy()
    df.reset_index(drop=True, inplace=True)

    
    
    if type == "fes":
        keep = ["3FE4S", "4FE4S", "4FE3S"]
        df = df[df.formula.isin(keep)]
        
        Y = df["formula"]

        if "formula" in df.columns:
            df = df.drop(columns="formula")
        if "id" in df.columns:
            df = df.drop(columns="id")
        if "smiles" in df.columns:
            df = df.drop(columns="smiles")
        if "res_name" in df.columns:
            df = df.drop(columns="res_name")
        if "type" in df.columns:
            df = df.drop(columns="type")
        if "aug_dist" in df.columns:
            df = df.drop(columns="aug_dist")
        if "structure" in df.columns:
            df = df.drop(columns="structure")
        
        X = df

    else:
        X = df.drop("formula", axis=1)
        X = X.fillna(0)
        y = df["formula"]
        Y = []
        if type == "as":
            for i in y:
                if i != "protein" and i != "active_site":
                    Y.append("protein")
                else:
                    Y.append(i)
        else:
            for i in y:
                if i == "active_site":
                    Y.append("protein")
                elif i == "3FE4S" or i == "4FE4S" or i == "4FE3S":
                    Y.append("fes_pocket")
                else:
                    Y.append("protein")  
        
        Y = pd.Series(Y)
    return X, Y

def train_model(
    directory: Path,
    finerprinttsv_name: str,
    model_name: str,
    model_type: str,
    model_dir_name: str,
    training_params: dict,
    fingerprint_type : dict,
    f_radius : float,
    pred_type : str
):
    """ """
    for subset in os.listdir(directory):
        print("train model in", subset, "...")
        subset_path = directory / Path(str(subset))
        fingeprint_db_path = subset_path / Path("db") / Path("fingerprintDB")
        if model_dir_name != None:
            model_dir_name_path = fingeprint_db_path / Path(model_dir_name)
            os.makedirs(model_dir_name_path, exist_ok=True)
        else:
            model_dir_name_path = fingeprint_db_path
        fingerprint_tsv_path = fingeprint_db_path / Path(finerprinttsv_name)

        
        X, Y = prepare_data(fingerprint_tsv_path, pred_type)
        model = build_model(model_type)

        if pred_type == "fes":
            fes_model = get_pipeline(model)
            fes_model.set_params(**training_params)

            fes_model.fit(X, Y)

            fes_model_output_path = model_dir_name_path / Path("fes_type_" + model_name + ".pkl")

            save_model(fes_model_output_path, fes_model, training_params, fingerprint_type, f_radius)
        elif pred_type == "fes_pocket":
            fes_pocket_model = get_pipeline(model)
            fes_pocket_model.set_params(**training_params)

            fes_pocket_model.fit(X, Y)

            fes_pocket_model_output_path = model_dir_name_path / Path("fes_pocket_" + model_name + ".pkl")

            save_model(fes_pocket_model_output_path, fes_pocket_model, training_params, fingerprint_type, f_radius)
        else:
            # define model
            as_model = get_pipeline(model)
            as_model.set_params(**training_params)

            # fit model
            as_model.fit(X, Y)
            
            # save model
            ## if modelname already exists override
            as_model_output_path = model_dir_name_path / Path("as_" + model_name + ".pkl")

            save_model(as_model_output_path, as_model, training_params, fingerprint_type, f_radius)

    if pred_type == "fes":
        return Path(model_dir_name)/Path("fes_type_" + model_name + ".pkl")
    elif pred_type == "fes_pocket":
        return Path(model_dir_name)/Path("fes_pocket_" + model_name + ".pkl")
    else:
        return Path(model_dir_name)/Path("as_" + model_name + ".pkl")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("dir", type=Path)
    parser.add_argument("fingerprinttsv_name", type=str)
    parser.add_argument("model_name", type=str)
    parser.add_argument("--model_dir_name", type=str, default=None)
    parser.add_argument(
        "--parameter_grid", action="store_true"
    )  # save all models with name_00X and save a extra file that contains name-param pairs
    parser.add_argument("--iter", type=int, default=1)
    parser.add_argument("--type", choices=["fes", "as", "fes_pocket"], default="fes")

    args = parser.parse_args()

    if args.parameter_grid:
        if args.model_dir_name == None:
            raise ValueError("Need to define model_dir_name")
        train_models(
            args.dir,
            args.fingerprinttsv_name,
            args.model_name,
            args.model_dir_name,
            args.iter,
        )
    else:
        train_model(
            args.dir,
            args.fingerprinttsv_name,
            args.model_name,
            args.model_dir_name,
            PARAMS,
            args.type
        )

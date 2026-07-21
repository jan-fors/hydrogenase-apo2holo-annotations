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

from src.utils.constants import COFACTOR_BLACKLIST
warnings.filterwarnings("ignore", module="sklearn")

from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression

RANDOM_STATE = 161
# CAUTION: if svm is used include: probablity = True
MODEL = MLPClassifier(random_state=RANDOM_STATE)  # insert model
PARAMS = {
    'clf__activation': 'relu',
    'clf__alpha': np.float64(0.0072349310782850096),
    'clf__batch_size': 128,
    'clf__beta_1': np.float64(0.9014878985254476),
    'clf__beta_2': np.float64(0.9911840289428955),
    'clf__hidden_layer_sizes': (256, 128, 64),
    'clf__learning_rate': 'invscaling',
    'clf__learning_rate_init': np.float64(0.008589062084956335),
    'clf__max_iter': 708,
    'clf__momentum': np.float64(0.7282085760409726),
    'clf__solver': 'adam',
    'clf__tol': np.float64(3.200280851029404e-05),
    'clf__validation_fraction': np.float64(0.1608894257116988)
    }

# PARAMS = {'clf__activation': 'tanh',                                                                                                                                
#  'clf__alpha': np.float64(0.0007131723391696173),                                                                                                          
#  'clf__batch_size': 32,                                                                                                                                    
#  'clf__beta_1': np.float64(0.9577537785150572),                                                                                                            
#  'clf__beta_2': np.float64(0.9336881675144799),                                                                                                            
#  'clf__hidden_layer_sizes': (1032, 256),                                                                                                                   
#  'clf__learning_rate': 'adaptive',                                                                                                                         
#  'clf__learning_rate_init': np.float64(0.002905802905049201),                                                                                              
#  'clf__max_iter': 773,                                                                                                                                     
#  'clf__momentum': np.float64(0.7996025744376337),                                                                                                          
#  'clf__solver': 'sgd',                                                                                                                                     
#  'clf__tol': np.float64(1.5278636601050295e-05),                                                                                                           
#  'clf__validation_fraction': np.float64(0.18966791991807014)}   

# MODEL = SVC(probability=True)
# PARAMS = {'clf__C': np.float64(14.011699774992735),
#  'clf__class_weight': None,
#  'clf__coef0': np.float64(6.5287809780266),
#  'clf__degree': 3,
#  'clf__gamma': np.float64(0.05158936159969933),
#  'clf__kernel': 'rbf',
#  'clf__shrinking': True,
#  'clf__tol': np.float64(0.0010111861961050502)}

# MODEL = RandomForestClassifier()
# PARAMS = {'clf__bootstrap': True,
#  'clf__class_weight': 'balanced',
#  'clf__max_depth': 20,
#  'clf__max_features': 'sqrt',
#  'clf__max_samples': np.float64(0.7511907372444331),
#  'clf__min_samples_leaf': 2,
#  'clf__min_samples_split': 3,
#  'clf__n_estimators': 116}

PARAM_DISTRIBUTION = {
    
}

# FUNCTIONS
def printl(text: str):
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{now}] {text}")

def get_pipeline():
    """ """
    # define model
    pipe = Pipeline(
        [
            ("scaler", StandardScaler()),
            ("clf", MODEL),
        ]
    )
    return pipe


def draw_params():
    sample = {}
    for key, dist in PARAM_DISTRIBUTION.items():
        if isinstance(dist, list):
            sample[key] = dist[np.random.randint(len(dist))]
        else:
            sample[key] = dist.rvs()
    return sample


def draw_unique(drawn, max_tries=100):
    for _ in range(max_tries):
        sample = draw_params()
        if sample not in drawn:
            drawn.append(sample)
            return sample
    raise ValueError("Could not draw a unique sample after max_tries")


def train_models(
    directory: Path,
    fingerprinttsv_name: str,
    model_name: str,
    model_dir_name: str,
    iter: int,
):
    previous_params = []
    for i in range(iter):
        print("Iteration", i)
        params = draw_unique(previous_params)

        i_model_name = model_name + "__" + str(i)

        train_model(
            directory, fingerprinttsv_name, i_model_name, model_dir_name, params
        )

def prepare_data(data: Path, random_state: int, type : str):
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
    model_dir_name: str,
    params: dict,
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

        #df = pd.read_csv(fingerprint_tsv_path, sep="\t")
        
        X, Y = prepare_data(fingerprint_tsv_path, RANDOM_STATE, pred_type)

        if pred_type == "fes":
            fes_model = get_pipeline()
            fes_model.set_params(**params)

            fes_model.fit(X, Y)

            fes_model_output_path = model_dir_name_path / Path("fes_type_" + model_name + ".pkl")
            with open(fes_model_output_path, "wb") as f:
                pickle.dump(fes_model, f)
        elif pred_type == "fes_pocket":
            fes_pocket_model = get_pipeline()
            fes_pocket_model.set_params(**params)

            fes_pocket_model.fit(X, Y)

            fes_pocket_model_output_path = model_dir_name_path / Path("fes_pocket_" + model_name + ".pkl")
            with open(fes_pocket_model_output_path, "wb") as f:
                pickle.dump(fes_pocket_model, f)
        else:
            # define model
            as_model = get_pipeline()
            as_model.set_params(**params)

            # fit model
            as_model.fit(X, Y)
            
            # save model
            ## if modelname already exists override
            as_model_output_path = model_dir_name_path / Path("as_" + model_name + ".pkl")
            with open(as_model_output_path, "wb") as f:
                pickle.dump(as_model, f)


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

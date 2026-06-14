"""
Creates a specific model for each training dataset

Create a list of params for the models and generate all of them to check afterwards

Allow to generate multiple models sampling from a distribution of parameters for RandomizedCV
"""
import os
from pathlib import Path
import argparse
import pandas as pd
import numpy as np
import pickle
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from scipy.stats import loguniform, uniform, randint

import warnings
warnings.filterwarnings("ignore", module="sklearn")

from sklearn.neural_network import MLPClassifier
from sklearn.ensemble import RandomForestClassifier

from sklearn.svm import SVC

RANDOM_STATE = 161
MODEL = SVC(random_state=RANDOM_STATE, probability=True) # insert model
PARAMS = {
    # fill with model specific params. Example: 
    # "clf__solver": "adam"
}

PARAM_DISTRIBUTION = {
    # "clf__hidden_layer_sizes":  [(64,), (128,), (256,),
    #                              (64, 32), (128, 64), (256, 128),
    #                              (128, 64, 32), (256, 128, 64),
    #                              (512, 256, 128)],
    # "clf__activation":          ["relu", "relu", "relu", "tanh", "logistic"],
    # "clf__solver":              ["adam", "adam", "adam", "sgd"],
    # "clf__alpha":               loguniform(1e-5, 1e-1),
    # "clf__learning_rate":       ["constant", "adaptive", "invscaling"],
    # "clf__learning_rate_init":  loguniform(1e-4, 1e-1),
    # "clf__max_iter":            randint(200, 1000),
    # "clf__tol":                 loguniform(1e-5, 1e-2),
    # "clf__early_stopping":      [False],
    # "clf__validation_fraction": uniform(0.1, 0.2),
    # "clf__batch_size":          [32, 64, 128, 256, "auto"],
    # "clf__momentum":            uniform(0.5, 0.45),
    # "clf__beta_1":              uniform(0.85, 0.14),
    # "clf__beta_2":              uniform(0.9, 0.099)
    # "clf__n_estimators": randint(100, 1000),
    #     "clf__max_depth": [None, 5, 10, 20, 30, 50],
    #     "clf__min_samples_split": randint(2, 20),
    #     "clf__min_samples_leaf": randint(1, 10),
    #     "clf__max_features": ["sqrt", "log2", None],
    #     "clf__class_weight": [None, "balanced"],
    #     "clf__bootstrap": [True, False],
    #   #  "clf__max_samples": uniform(0.5, 0.5),
      "clf__C":               loguniform(1e-3, 1e4),
    "clf__kernel":          ["rbf", "rbf", "rbf", "linear", "poly", "sigmoid"],
    "clf__gamma":           ["scale", "auto"] + list(loguniform(1e-4, 1e1).rvs(20)),
    "clf__degree":          randint(2, 6),
    "clf__coef0":           uniform(0, 10),
    "clf__class_weight":    [None, "balanced"],
    "clf__tol":             loguniform(1e-5, 1e-2),
    "clf__shrinking":       [True, False],
}

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

def train_models(directory : Path, fingerprinttsv_name : str, model_name : str, model_dir_name : str, iter : int):
    previous_params = []
    for i in range(iter):
        print("Iteration", i)
        params = draw_unique(previous_params)

        i_model_name = model_name + "__" + str(i)

        train_model(directory, fingerprinttsv_name, i_model_name, model_dir_name, params)

def train_model(directory : Path, finerprinttsv_name : str, model_name : str, model_dir_name : str,  params : dict):
    for subset in os.listdir(directory):
        print("train model in", subset, "...")
        subset_path = directory / Path(str(subset))
        fingeprint_db_path = subset_path / Path('db') / Path('fingerprintDB')
        if model_dir_name != None:
            model_dir_name_path = fingeprint_db_path / Path(model_dir_name)
            os.makedirs(model_dir_name_path, exist_ok=True)
        else:
            model_dir_name_path = fingeprint_db_path
        fingerprint_tsv_path = fingeprint_db_path / Path(finerprinttsv_name)
        
        
        db = pd.read_csv(fingerprint_tsv_path, sep="\t")
        # prep X and y
        X = db.drop(columns=['formula', 'id', 'smiles', 'res_name'])
        X_np = X.to_numpy()

        Y = db["formula"]
        
        # define model
        model = get_pipeline()
        model.set_params(**params)

        # fit model
        model.fit(X, Y)

        # save model 
        ## if modelname already exists override
        model_output_path = model_dir_name_path / Path(model_name + ".pkl")
        with open(model_output_path, "wb") as f:
            pickle.dump(model, f)



if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument("dir", type=Path)
    parser.add_argument("fingerprinttsv_name", type=str)
    parser.add_argument("model_name", type=str)
    parser.add_argument("--model_dir_name", type=str, default=None)
    parser.add_argument("--parameter_grid", action="store_true") # save all models with name_00X and save a extra file that contains name-param pairs
    parser.add_argument("--iter", type=int, default=1)

    args = parser.parse_args()

    if args.parameter_grid:
        if args.model_dir_name == None:
            raise ValueError("Need to define model_dir_name")
        train_models(args.dir, args.fingerprinttsv_name, args.model_name, args.model_dir_name, args.iter)
    else:
        train_model(args.dir, args.fingerprinttsv_name, args.model_name, args.model_dir_name, PARAMS)

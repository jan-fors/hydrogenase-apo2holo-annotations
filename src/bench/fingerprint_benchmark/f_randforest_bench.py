""" """

# IMPORTS
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.experimental import enable_halving_search_cv
from sklearn.model_selection import RandomizedSearchCV
from pathlib import Path
import pandas as pd
from pprint import pprint
from src.utils.constants import COFACTOR_BLACKLIST
import datetime
from sklearn.ensemble import RandomForestClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import confusion_matrix, classification_report
from scipy.stats import loguniform, randint, uniform
import numpy as np
import warnings
import argparse
import time
import pickle

warnings.filterwarnings("ignore", module="sklearn")


# FUNCTIONS
def printl(text: str):
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{now}] {text}")


def get_pipeline(random_state: int):
    """ """
    # define model
    pipe = Pipeline(
        [
            ("scaler", StandardScaler()),
            ("clf", RandomForestClassifier(random_state=random_state)),
        ]
    )
    return pipe


def get_param_distributions():
    param_distributions_rf = {
        "clf__n_estimators": randint(100, 1000),
        "clf__max_depth": [None, 5, 10, 20, 30, 50],
        "clf__min_samples_split": randint(2, 20),
        "clf__min_samples_leaf": randint(1, 10),
        "clf__max_features": ["sqrt", "log2", None],
        "clf__class_weight": [None, "balanced"],
        "clf__bootstrap": [True, False],
        "clf__max_samples": uniform(0.5, 0.5),
    }
    return param_distributions_rf


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


def bench(data, random_state, test_size, scoring, jobs, k, n_iter, type, save_model) -> dict:
    """ """
    result = {}

    X, y = prepare_data(data, random_state, type)

    # initial split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state, stratify=y
    )

    # train
    search = RandomizedSearchCV(
        get_pipeline(random_state),
        get_param_distributions(),
        n_iter=n_iter,
        scoring=scoring,
        cv=k,
        random_state=random_state,
        n_jobs=jobs,
    ).fit(X_train, y_train)

    result["best_params"] = search.best_params_

    # k-fold cv on training data
    scores = cross_val_score(
        search.best_estimator_,  # pipeline with best params already baked in
        X_train,
        y_train,
        cv=k,
        scoring=scoring,
        n_jobs=jobs,
    )

    result["cv_mean"] = scores.mean()
    result["cv_std"] = scores.std()

    # final test on testset
    final_model = search.best_estimator_
    y_pred = final_model.predict(X_test)

    result["classification_report"] = classification_report(y_test, y_pred, digits=6)

    if save_model:
        # save the final model + metadata
        dataname = data.stem
        out_dir = data.parent / Path("models")
        out_dir.mkdir(exist_ok=True)
        stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        model_path = out_dir / f"mlp_{type}_data_{dataname}_rs{random_state}_{stamp}.pkl"
        with open(model_path, "wb") as f:
            pickle.dump(final_model, f)
        printl(f"Saved model + metadata to {model_path}")

    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("data",type=Path,  help="Path to the fingerprint tsv file.")
    parser.add_argument(
        "--random-state", default=161, type=int, help="Random state for reproducibility"
    )
    parser.add_argument(
        "--test-size",
        default=0.1,
        type=float,
        help="Relative size of initially taken test samples",
    )
    parser.add_argument(
        "--scoring", choices=["accuracy", "f1_weighted", "f1_macro"], default="f1_macro"
    )
    parser.add_argument(
        "--jobs", type=int, default=1, help="Number of parallel Threads. [1]"
    )
    parser.add_argument(
        "--k", type=int, default=5, help="Amount of Cross Validation Rounds. [5]"
    )
    parser.add_argument("--n-iter", type=int, default=10, help="Number of iterations. [10]")
    parser.add_argument("--type", choices=["fes", "as", "fes_pocket"], default="as")
    args = parser.parse_args()

    bench(
        args.data,
        args.random_state,
        args.test_size,
        args.scoring,
        args.jobs,
        args.k,
        args.n_iter,
        args.type
    )

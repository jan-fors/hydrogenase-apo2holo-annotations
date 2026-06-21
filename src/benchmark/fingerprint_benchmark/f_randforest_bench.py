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


def prepare_data(data: Path, test_size: float, random_state: int):
    """ """
    df = pd.read_csv(data, sep="\t")
    printl("Read data...")

    # df.drop("id", axis=1, inplace=True)
    # df.drop("res_name", axis=1, inplace=True)
    # df.drop("smiles", axis=1, inplace=True)
    df = df[~df["formula"].isin(COFACTOR_BLACKLIST)].copy()
    df.reset_index(drop=True, inplace=True)

    X = df.drop("formula", axis=1)
    y = df["formula"]

    return X, y


def bench(data, random_state, test_size, scoring, jobs, k, n_iter):
    """ """
    printl(f"Random state={random_state}")
    printl(f"Test Size={test_size}")
    printl(f"Scoring={scoring}")
    printl(f"CV={k}")
    printl(f"Iterations={n_iter}")

    X, y = prepare_data(data, test_size, random_state)

    # initial split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state, stratify=y
    )

    printl(
        f"Create initial train/test-split with test-size={test_size} and random_state={random_state}"
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

    printl(f"Best params are:")
    pprint(search.best_params_)
    print("\n")

    # k-fold cv on training data
    scores = cross_val_score(
        search.best_estimator_,  # pipeline with best params already baked in
        X_train,
        y_train,
        cv=k,
        scoring=scoring,
    )

    printl(f"CV {scoring}: {scores.mean():.3f} (+/- {scores.std():.3f})")

    # final test on testset
    final_model = search.best_estimator_
    y_pred = final_model.predict(X_test)
    print(classification_report(y_test, y_pred, digits=3))

    # confusion matrix
    cm = confusion_matrix(y_test, y_pred)
    labels = sorted(y_test.unique())

    pad = max(len(str(l)) for l in labels) + 2

    print("\nConfusion Matrix:")
    print(f"{'':>{pad}}", end="")
    for label in labels:
        print(f"{str(label):>{pad}}", end="")
    print()
    for i, row_label in enumerate(labels):
        print(f"{str(row_label):>{pad}}", end="")
        for val in cm[i]:
            print(f"{val:>{pad}}", end="")
        print()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("data", help="Path to the fingerprint tsv file.")
    parser.add_argument(
        "--random-state", default=0, type=int, help="Random state for reproducibility"
    )
    parser.add_argument(
        "--test-size",
        default=0.1,
        type=float,
        help="Relative size of initially taken test samples",
    )
    parser.add_argument(
        "--scoring", choices=["accuracy", "f1_weighted", "f1_macro"], default="accuracy"
    )
    parser.add_argument(
        "--jobs", type=int, default=1, help="Number of parallel Threads"
    )
    parser.add_argument(
        "--k", type=int, default=5, help="Amount of Cross Validation Rounds"
    )
    parser.add_argument("--n-iter", type=int, default=10, help="Number of iterations")

    args = parser.parse_args()

    bench(
        args.data,
        args.random_state,
        args.test_size,
        args.scoring,
        args.jobs,
        args.k,
        args.n_iter,
    )

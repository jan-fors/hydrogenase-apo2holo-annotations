"""
Script that performs a HalvinRandomSearchCV to find the best params and check with a k-fold cv on the training dataset.
Returns a final accuracy on an unseen testset
"""

# IMPORTS
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.experimental import enable_halving_search_cv
from sklearn.model_selection import RandomizedSearchCV
from pathlib import Path
import pandas as pd
import sys
from pprint import pprint
sys.path.append("..")
from src.utils.constants import COFACTOR_BLACKLIST
import datetime
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import confusion_matrix, classification_report
import argparse
from sklearn.neural_network import MLPClassifier
from scipy.stats import loguniform, uniform, randint
import warnings

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
        ("clf", MLPClassifier(random_state=random_state)),  
    ]
    )
    return pipe


def get_param_distributions():
    param_distributions_mlp = {
    "clf__hidden_layer_sizes":  [(64,), (128,), (256,),
                                 (64, 32), (128, 64), (256, 128),
                                 (128, 64, 32), (256, 128, 64),
                                 (512, 256, 128)],
    "clf__activation":          ["relu", "relu", "relu", "tanh", "logistic"],
    "clf__solver":              ["adam", "adam", "adam", "sgd"],
    "clf__alpha":               loguniform(1e-5, 1e-1),
    "clf__learning_rate":       ["constant", "adaptive", "invscaling"],
    "clf__learning_rate_init":  loguniform(1e-4, 1e-1),
    "clf__max_iter":            randint(200, 1000),
    "clf__tol":                 loguniform(1e-5, 1e-2),
    "clf__early_stopping":      [False],
    "clf__validation_fraction": uniform(0.1, 0.2),
    "clf__batch_size":          [32, 64, 128, 256, "auto"],
    "clf__momentum":            uniform(0.5, 0.45),
    "clf__beta_1":              uniform(0.85, 0.14),
    "clf__beta_2":              uniform(0.9, 0.099)
    }
    return param_distributions_mlp


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

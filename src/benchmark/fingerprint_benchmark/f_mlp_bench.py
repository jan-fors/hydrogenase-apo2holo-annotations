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
        ("clf", MLPClassifier(random_state=random_state)),  
    ]
    )
    return pipe


def get_param_distributions():
    param_distributions_mlp = {
    "clf__hidden_layer_sizes": [
        # single layer, small to large
        (20,), (40,), (64,), (128,), (256,), (512,), (1032,), (2064,),
        # two-layer funnels
        (40, 20), (64, 32), (128, 64), (256, 128), (512, 256), (1032, 516),
        (2064, 512), (1032, 256),
        # three-layer funnels
        (128, 64, 32), (256, 128, 64), (512, 256, 128), (2064, 1032, 516),
        (512, 128, 32),
        # four-layer funnels
        (512, 256, 128, 64), (2064, 1032, 516, 128),
        # symmetric / expand-contract
        (20, 40, 20), (64, 128, 64), (128, 256, 128),
        (20, 40, 60, 40, 20), (64, 128, 256, 128, 64),
        # aggressive reduction for high-dim input (5000 → small)
        (2064, 512, 128), (2064, 256), (1032, 128, 16),
    ],
      "clf__activation":          ["relu", "tanh", "logistic"],  # weight via separate means if needed
    "clf__solver":              ["adam", "sgd"],
    "clf__alpha":               loguniform(1e-5, 1e-1),
    "clf__learning_rate":       ["constant", "adaptive", "invscaling"],
    "clf__learning_rate_init":  loguniform(1e-4, 1e-2),
    "clf__max_iter":            randint(200, 1000),
    "clf__tol":                 loguniform(1e-5, 1e-2),
#    "clf__early_stopping":      [True],
    "clf__validation_fraction": uniform(0.1, 0.2),
    "clf__batch_size":          [32, 64, 128, 256, "auto"],
    "clf__momentum":            uniform(0.5, 0.45),
    "clf__beta_1":              uniform(0.85, 0.14),
    "clf__beta_2":              uniform(0.9, 0.099),
    }
    return param_distributions_mlp


def prepare_data(data: Path, test_size: float, random_state: int, type : str):
    """ """
    df = pd.read_csv(data, sep="\t")
    printl("Read data...")

    if "id" in df.columns:
        df.drop("id", axis=1, inplace=True)
    if "res_name" in df.columns:
        df.drop("res_name", axis=1, inplace=True)
    if "smiles" in df.columns:
        df.drop("smiles", axis=1, inplace=True)

    if "fingerprint" in df.columns:
        df.drop("fingerprint", axis=1, inplace=True)
    df = df[~df["formula"].isin(COFACTOR_BLACKLIST)].copy()
    df.reset_index(drop=True, inplace=True)

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
            else:
                Y.append(i)  
    
    Y = pd.Series(Y)
    return X, Y


def bench(data, random_state, test_size, scoring, jobs, k, n_iter, type):
    """ """
    t_start = time.perf_counter()

    printl(f"Random state={random_state}")
    printl(f"Test Size={test_size}")
    printl(f"Scoring={scoring}")
    printl(f"CV={k}")
    printl(f"Iterations={n_iter}")

    X, y = prepare_data(data, test_size, random_state, type)

    # initial split
    X_train, X_test, y_train, y_test = train_test_split(
        X, y, test_size=test_size, random_state=random_state, stratify=y
    )

    printl(
        f"Create initial train/test-split with test-size={test_size} and random_state={random_state}"
    )

    # train
    t0 = time.perf_counter()
    search = RandomizedSearchCV(
        get_pipeline(random_state),
        get_param_distributions(),
        n_iter=n_iter,
        scoring=scoring,
        cv=k,
        random_state=random_state,
        n_jobs=jobs,
    ).fit(X_train, y_train)
    t_search = time.perf_counter() - t0
    printl(f"RandomizedSearchCV took {t_search:.1f}s ({t_search/60:.1f} min)")

    printl(f"Best params are:")
    pprint(search.best_params_)
    print("\n")

    # k-fold cv on training data
    t0 = time.perf_counter()
    scores = cross_val_score(
        search.best_estimator_,  # pipeline with best params already baked in
        X_train,
        y_train,
        cv=k,
        scoring=scoring,
        n_jobs=jobs,
    )
    t_cv = time.perf_counter() - t0
    printl(f"cross_val_score took {t_cv:.1f}s ({t_cv/60:.1f} min)")

    printl(f"CV {scoring}: {scores.mean():.3f} (+/- {scores.std():.3f})")

    # final test on testset
    t0 = time.perf_counter()
    final_model = search.best_estimator_
    y_pred = final_model.predict(X_test)
    t_predict = time.perf_counter() - t0
    printl(f"Final predict took {t_predict:.2f}s")

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

    t_total = time.perf_counter() - t_start
    printl(f"TOTAL runtime: {t_total:.1f}s ({t_total/60:.1f} min)")

    # save the final model + metadata
    dataname = data.stem
    out_dir = data.parent / Path("models")
    out_dir.mkdir(exist_ok=True)
    stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
    model_path = out_dir / f"mlp_{type}_data{dataname}_rs{random_state}_{stamp}.pkl"

    with open(model_path, "wb") as f:
        pickle.dump(final_model, f)
    printl(f"Saved model + metadata to {model_path}")

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
        "--jobs", type=int, default=1, help="Number of parallel Threads"
    )
    parser.add_argument(
        "--k", type=int, default=5, help="Amount of Cross Validation Rounds"
    )
    parser.add_argument("--n-iter", type=int, default=10, help="Number of iterations")
    parser.add_argument("--type", choices=["fes", "as"], default="as")
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

from pathlib import Path
from typing import Literal

from apo2holo.bench.fingerprint_benchmark.f_logreg_bench import bench as lr_bench
from apo2holo.bench.fingerprint_benchmark.f_mlp_bench import bench as mlp_bench
from apo2holo.bench.fingerprint_benchmark.f_randforest_bench import bench as rf_bench
from apo2holo.bench.fingerprint_benchmark.f_svm_bench import bench as svm_bench
from apo2holo.bench.fingerprint_benchmark.f_knn_bench import bench as knn_bench
from apo2holo.bench.fingerprint_benchmark.f_lda_bench import bench as lda_bench
from apo2holo.bench.fingerprint_benchmark.f_qda_bench import bench as qda_bench

MODEL_TYPE_REGISTRY = {
    "mlp":mlp_bench,
    "rf":rf_bench,
    "svm":svm_bench,
    "lr":lr_bench,
    "knn": knn_bench, # https://scikit-learn.org/stable/modules/generated/sklearn.neighbors.KNeighborsClassifier.html
    "lda": lda_bench, # https://scikit-learn.org/stable/modules/generated/sklearn.discriminant_analysis.LinearDiscriminantAnalysis.html 
    "qda": qda_bench # https://scikit-learn.org/stable/modules/generated/sklearn.discriminant_analysis.QuadraticDiscriminantAnalysis.html
}

def perform_gridsearch(
    model_type : Literal["mlp", "rf", "svm", "lr", "knn", "lda", "qda"],
    data: Path,
    random_state: int = 161,
    test_size: float = 0.1,
    scoring: Literal["accuracy", "f1_weighted", "f1_macro"] = "f1_macro",
    jobs: int = 1,
    k: int = 5,
    n_iter: int = 10,
    type=["fes", "as", "fes_pocket"],
    save_model: bool = False,
    do_pca : bool = False
):
    """
    """
    model_bench = MODEL_TYPE_REGISTRY[model_type]

    return model_bench(data, random_state, test_size, scoring, jobs, k, n_iter, type, save_model, do_pca)
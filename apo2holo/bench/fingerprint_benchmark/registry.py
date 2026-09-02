from pathlib import Path
from typing import Literal

from apo2holo.bench.fingerprint_benchmark.f_logreg_bench import bench as lr_bench
from apo2holo.bench.fingerprint_benchmark.f_mlp_bench import bench as mlp_bench
from apo2holo.bench.fingerprint_benchmark.f_randforest_bench import bench as rf_bench
from apo2holo.bench.fingerprint_benchmark.f_svm_bench import bench as svm_bench

MODEL_TYPE_REGISTRY = {
    "mlp":mlp_bench,
    "rf":rf_bench,
    "svm":svm_bench,
    "lr":lr_bench,
}

def perform_gridsearch(
    model_type : Literal["mlp", "rf", "svm", "lr"],
    data: Path,
    random_state: int = 161,
    test_size: float = 0.1,
    scoring: Literal["accuracy", "f1_weighted", "f1_macro"] = "f1_macro",
    jobs: int = 1,
    k: int = 5,
    n_iter: int = 10,
    type=["fes", "as", "fes_pocket"],
    save_model: bool = False,
):
    """
    """
    model_bench = MODEL_TYPE_REGISTRY[model_type]

    return model_bench(data, random_state, test_size, scoring, jobs, k, n_iter, type, save_model)
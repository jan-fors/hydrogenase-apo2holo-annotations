import pickle
from pathlib import Path

def load_model(model_path : Path):
    """
    Returns model : model obj, training_params : dict, fingerprint_type : dict, f_radius : float
    """

    with open(model_path, "rb") as f:
        data = pickle.load(f)

    return data["model"], data["training_params"], data["fingerprint_type"], data["f_radius"]
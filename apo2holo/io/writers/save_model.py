import pickle
from pathlib import Path

def save_model(output_path : Path, model, training_params, fingerprint_type, f_radius):
    """
    Part of the model output dir should be: 
    1. the model itself
    2. the training params
    3. the corresponding fingerprints
    4. the f_radius
    """
    data = {
        "model": model,
        "training_params": training_params,
        "fingerprint_type": fingerprint_type,
        "f_radius": f_radius
    }
    with open(output_path, "wb") as f:
        pickle.dump(data, f)
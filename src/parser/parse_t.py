import numpy as np

def parse_t(t_str: str) -> np.ndarray:
    vals = [float(x) for x in t_str.split(",")]
    if len(vals) != 3:
        raise ValueError("t must have 3 comma-separated floats")
    return np.array(vals, dtype=float)
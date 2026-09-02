import numpy as np


def parse_u(u_str: str) -> np.ndarray:
    vals = [float(x) for x in u_str.split(",")]
    if len(vals) != 9:
        raise ValueError("u must have 9 comma-separated floats")
    return np.array(vals, dtype=float).reshape(3, 3)
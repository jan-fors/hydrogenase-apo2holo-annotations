import numpy as np

def calculate_geometric_centers(cofactors : dict[str, list[tuple[str, float, float, float]]]) -> dict[str, tuple[float, float, float]]:
    """
    """
    res = {}
    for key in cofactors.keys():
        atoms = cofactors[key]
        P = np.array([[a[1], a[2], a[3]] for a in atoms], dtype=float)
        res[key] = _calculate_center(P)

    return res


def _calculate_center(M : np.ndarray) -> tuple[float,float,float]:
    """
    """
    if M.ndim != 2 or M.shape[1] != 3:
        raise ValueError(f"Expected array of shape (N,3), got {M.shape}")

    center = np.mean(M, axis=0)

    return float(center[0]), float(center[1]), float(center[2])

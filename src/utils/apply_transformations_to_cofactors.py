import numpy as np

def apply_transformations_to_cofactors(cofactors : dict[str, list[tuple[str, float, float, float]]], u : np.array, t : np.array) -> dict[str, list[tuple[str, float, float, float]]]:
    """"""
    R = np.asarray(u, dtype=float).reshape(3, 3)
    for key in cofactors.keys():
        atoms = cofactors[key]
        P = np.array([[a[1], a[2], a[3]] for a in atoms], dtype=float).T
        out_array = _apply_rt_cols(P, R, t)
        out_array = out_array.T

        for i in range(len(cofactors[key])):
            cofactors[key][i] = (cofactors[key][i][0], out_array[i][0], out_array[i][1], out_array[i][2])
    
    return cofactors


def _apply_rt_cols(points_3xn: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    """
    points_3xn: (3,N)
    R: (3,3)
    t: (3,)
    returns: (3,N)
    """
    return R @ points_3xn + t[:, None]
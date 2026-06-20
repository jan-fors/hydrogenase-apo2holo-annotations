import numpy as np

def random_oversample(X, y, target_count=None, random_state=None, undersample=False):
    """
    Balance classes by randomly duplicating samples from minority classes.

    Parameters
    ----------
    X : array-like, shape (n_samples, n_features)
    y : array-like, shape (n_samples,)
    target_count : int or None
        Samples per class after resampling. If None, uses the size of the
        largest class (full balance).
    random_state : int or None

    Returns
    -------
    X_res, y_res : resampled arrays, shuffled.
    """
    rng = np.random.default_rng(random_state)
    X = np.asarray(X)
    y = np.asarray(y)

    classes, counts = np.unique(y, return_counts=True)
    if target_count is None:
        target_count = counts.max()

    idx_resampled = []
    for cls in classes:
        cls_idx = np.where(y == cls)[0]
        n = len(cls_idx)
        if n >= target_count:
            # keep all (or subsample if you set target below its size)
            if undersample:
                chosen = rng.choice(cls_idx, size=target_count, replace=False)
            else:
                chosen = cls_idx
            
        else:
            # keep originals, then sample the remainder WITH replacement
            extra = rng.choice(cls_idx, size=target_count - n, replace=True)
            chosen = np.concatenate([cls_idx, extra])
        idx_resampled.append(chosen)

    idx_resampled = np.concatenate(idx_resampled)
    rng.shuffle(idx_resampled)
    return X[idx_resampled], y[idx_resampled]

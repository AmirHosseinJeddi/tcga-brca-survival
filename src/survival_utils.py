"""Small helpers shared by the notebooks: fast univariate Cox screening and fixed CV folds."""
import json
import numpy as np
from scipy.stats import chi2
from sklearn.model_selection import RepeatedStratifiedKFold


def cox_score_screen(X, time, event):
    """Univariate Cox score test for every column of X at once (fast).

    For each gene, tests whether expression is associated with the hazard,
    using the Cox partial-likelihood score statistic at beta = 0 (Breslow ties).
    Returns (z, p). The sign of z shows direction: z > 0 means higher expression
    goes with higher risk of death. Rank genes by abs(z).
    """
    X = np.asarray(X, dtype=np.float64)
    time = np.asarray(time, dtype=np.float64)
    event = np.asarray(event).astype(bool)

    order = np.argsort(-time, kind="stable")          # latest time first
    t, e, Xs = time[order], event[order], X[order]
    n = len(t)

    # for each subject: position of the last member of its tie group
    # (risk set = everyone with time >= t_i = rows 0..last)
    change = t[1:] != t[:-1]
    group_id = np.concatenate([[0], np.cumsum(change)])
    group_end = np.concatenate([np.flatnonzero(change), [n - 1]])
    last = group_end[group_id]

    n_risk = (last + 1).astype(np.float64)
    S1 = np.cumsum(Xs, axis=0)[last]
    S2 = np.cumsum(Xs ** 2, axis=0)[last]

    ev = np.flatnonzero(e)
    mean = S1[ev] / n_risk[ev, None]
    var = S2[ev] / n_risk[ev, None] - mean ** 2

    U = (Xs[ev] - mean).sum(axis=0)
    V = var.sum(axis=0)
    ok = V > 1e-12
    z = np.zeros(X.shape[1])
    z[ok] = U[ok] / np.sqrt(V[ok])
    p = np.ones(X.shape[1])
    p[ok] = chi2.sf(z[ok] ** 2, 1)
    return z, p


def make_folds(strata, n_splits=5, n_repeats=10, seed=42):
    """Repeated stratified K-fold. Returns a list of (train_idx, test_idx)."""
    rskf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats, random_state=seed)
    return [(tr, te) for tr, te in rskf.split(np.zeros(len(strata)), strata)]


def save_folds(path, patients, folds):
    with open(path, "w") as f:
        json.dump({"patients": list(patients),
                   "folds": [{"train": tr.tolist(), "test": te.tolist()} for tr, te in folds]}, f)


def load_folds(path):
    with open(path) as f:
        d = json.load(f)
    folds = [(np.array(x["train"]), np.array(x["test"])) for x in d["folds"]]
    return d["patients"], folds

"""Exact check for the 37-cell core-loading problem.

Finds the proven minimum number of same-type adjacent pairs for a given
fuel inventory, using an exact MILP solver (HiGHS via SciPy).
Use it to certify the OpenJij energies.

Run:  python exact_check.py
"""
import time
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds

# Same 1-indexed edge list as In_core.ipynb
EDGES_1 = [
    (1,2),(1,5),(2,3),(2,6),(3,7),(4,5),(4,10),(5,6),(5,11),(6,7),
    (6,12),(7,8),(7,13),(8,14),(9,10),(9,16),(10,11),(10,17),(11,12),(11,18),
    (12,13),(12,19),(13,14),(13,20),(14,15),(14,21),(15,22),(16,17),(16,23),(17,18),
    (17,24),(18,19),(18,25),(19,20),(19,26),(20,21),(20,27),(21,22),(21,28),(22,29),
    (23,24),(24,25),(24,30),(25,26),(25,31),(26,27),(26,32),(27,28),(27,33),(28,29),
    (28,34),(30,31),(31,32),(31,35),(32,33),(32,36),(33,34),(33,37),(35,36),(36,37),
]
EDGES = [(u - 1, v - 1) for u, v in EDGES_1]
N = 37


def min_same_type(n_burnt: int):
    """Exact minimum of same-type adjacencies with n_burnt burnt assemblies.

    Variables: x_i = 1 if cell i is burnt; y_e = 1 if edge e is burnt-burnt;
    z_e = 1 if edge e is fresh-fresh. Standard linearization of x_u * x_v.
    """
    m = len(EDGES)
    nv = N + 2 * m
    c = np.r_[np.zeros(N), np.ones(2 * m)]
    rows, lb, ub = [], [], []
    for e, (u, v) in enumerate(EDGES):
        r = np.zeros(nv); r[u] = r[v] = 1; r[N + e] = -1          # x_u + x_v - 1 <= y_e
        rows.append(r); lb.append(-np.inf); ub.append(1)
        r = np.zeros(nv); r[u] = r[v] = -1; r[N + m + e] = -1     # 1 - x_u - x_v <= z_e
        rows.append(r); lb.append(-np.inf); ub.append(-1)
    r = np.zeros(nv); r[:N] = 1                                   # inventory
    rows.append(r); lb.append(n_burnt); ub.append(n_burnt)

    t = time.perf_counter()
    res = milp(c, constraints=LinearConstraint(np.array(rows), lb, ub),
               integrality=np.r_[np.ones(N), np.zeros(2 * m)],
               bounds=Bounds(0, 1))
    elapsed = time.perf_counter() - t
    layout = np.round(res.x[:N]).astype(int)
    return round(res.fun), elapsed, layout


if __name__ == "__main__":
    for n_burnt in [15, 16]:          # the two cases in the RET report
        best, secs, _ = min_same_type(n_burnt)
        print(f"{n_burnt:2d} burnt / {N - n_burnt:2d} fresh -> "
              f"proven minimum = {best} same-type adjacencies ({secs:.3f} s)")

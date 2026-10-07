"""MK8000 vs DWM1001C - UWB indoor positioning.

Position estimation pipeline of the paper, applied to the raw measurements in ``data/``.

    python uwb.py        # recompute Tables II and III from the raw data and check them against the paper

Pipeline, per test point (30 repeated two-way-ranging cycles per anchor):
  1. project every slant range onto the horizontal plane;
  2. average the 30 projected ranges per anchor;
  3. pick a reference anchor and, if more than four anchors are available, drop the weakest one;
  4. trilaterate with OLS, per-point WLS or RANSAC (linearised TOA system, Caffery).
"""
from __future__ import annotations

import itertools
import sys
from functools import lru_cache
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA, RESULTS = ROOT / "data", ROOT / "results"

# Table I - anchor coordinates (x, y, z) in metres; tag height Z_TAG.
ANCHORS = {
    "A0": (1.960, 0.400, 1.73),
    "A1": (4.965, 0.400, 1.73),
    "A2": (1.600, 11.737, 0.65),
    "A3": (5.965, 17.420, 0.65),
    "A4": (11.031, 15.995, 0.65),
}
Z_TAG = 1.00
# 64 interior test points + 4 exterior control points (IDs 1, 26, 112, 129), 30 cycles each.
FILES = {"DWM1001C": ("dwm1001_measurements.csv", "dwm1001_exterior.csv"),
         "MK8000": ("mk8000_measurements.csv", "mk8000_exterior.csv")}
INTERIOR = "punct_interior"
ESTIMATORS = ("OLS", "WLS", "RANSAC")

# Table III of the paper: (mean, p95) over the 64 interior points and mean over all 68 points, in cm.
PAPER_TABLE_III = {
    ("DWM1001C", "OLS"): (105.0, 254.1, 104.1), ("DWM1001C", "WLS"): (69.4, 148.7, 70.2),
    ("DWM1001C", "RANSAC"): (61.2, 137.3, 63.0),
    ("MK8000", "OLS"): (247.7, 587.0, 249.0), ("MK8000", "WLS"): (159.7, 467.9, 164.6),
    ("MK8000", "RANSAC"): (203.7, 681.8, 200.1),
}
# Table II of the paper: (a0, a1, eps, n) per system and anchor.
PAPER_TABLE_II = {
    ("DWM1001C", "A0"): (0.23, 0.99, 0.24, 49), ("DWM1001C", "A1"): (-0.16, 1.04, 0.60, 56),
    ("DWM1001C", "A2"): (-0.12, 1.02, 0.32, 68), ("DWM1001C", "A3"): (0.28, 1.04, 0.38, 68),
    ("DWM1001C", "A4"): (-0.14, 1.20, 1.57, 24),
    ("MK8000", "A0"): (-0.64, 1.17, 1.03, 68), ("MK8000", "A1"): (-1.35, 1.21, 1.28, 68),
    ("MK8000", "A2"): (-0.37, 1.06, 0.44, 68), ("MK8000", "A3"): (0.00, 1.04, 0.41, 68),
    ("MK8000", "A4"): (-0.40, 1.31, 1.48, 68),
}


# --------------------------------------------------------------------------- data
def load(system: str) -> pd.DataFrame:
    """All raw measurements of one system, with horizontally projected ranges ``h_A0..h_A4`` [m] added.
    Point 1 of the DWM1001C was measured six times; all of its 180 cycles are pooled."""
    df = pd.concat([pd.read_csv(DATA / f) for f in FILES[system]], ignore_index=True)
    for i, (a, (_, _, z)) in enumerate(ANCHORS.items()):
        d2 = df[f"r_AN{i}"] ** 2 - (z - Z_TAG) ** 2
        df[f"h_{a}"] = np.sqrt(d2.where(d2 > 0))          # NaN if missing or shorter than the height offset
    return df


# --------------------------------------------------------------------------- estimators
def _linearize(active, rng, ref):
    """A x = b from subtracting the equation of the reference anchor (x = tag position)."""
    rx, ry, _ = ANCHORS[ref]
    rows = [a for a in active if a != ref]
    A = np.array([[2 * (ANCHORS[a][0] - rx), 2 * (ANCHORS[a][1] - ry)] for a in rows])
    b = np.array([rng[ref] ** 2 - rng[a] ** 2 + ANCHORS[a][0] ** 2 - rx ** 2 + ANCHORS[a][1] ** 2 - ry ** 2
                  for a in rows])
    return A, b, rows


def ols(active, rng, ref):
    A, b, _ = _linearize(active, rng, ref)
    return np.linalg.lstsq(A, b, rcond=None)[0]


def wls(active, rng, ref, var, min_var=1e-6):
    """Weighted LS with W = diag(1/var_i); ``var`` = per-point variance of the repeated ranges [m^2]."""
    A, b, rows = _linearize(active, rng, ref)
    W = np.diag(1.0 / np.maximum([var[a] for a in rows], min_var))
    return np.linalg.solve(A.T @ W @ A, A.T @ W @ b)


def ransac(active, rng, tau):
    """Deterministic RANSAC: every 3-anchor subset is solved by OLS, the one with the largest consensus
    (| ||x - a_i|| - d_i | <= tau) wins, and the final position is an OLS fit on its inliers."""
    best, best_x = [], None
    for trio in itertools.combinations(active, 3):
        x = ols(trio, rng, trio[0])
        inliers = [a for a in active if abs(np.hypot(*(x - ANCHORS[a][:2])) - rng[a]) <= tau]
        if len(inliers) > len(best):
            best, best_x = inliers, x
    if len(best) < 3:
        return ols(active, rng, active[0]), list(active)
    return ols(best, rng, best[0]), best


def select_anchors(rng, quality, max_anchors=4):
    """Reference anchor = best quality (tie: the farther one). If more than ``max_anchors`` anchors are
    available, the weakest one is dropped (tie: the nearer one)."""
    avail = list(rng)
    ref = max(avail, key=lambda a: (quality[a], rng[a]))
    active = avail
    if len(avail) > max_anchors:
        weakest = min(avail, key=lambda a: (quality[a], rng[a]))
        active = [a for a in avail if a != weakest]
    return (ref if ref in active else active[0]), active


def quality(system: str, g: pd.DataFrame) -> dict:
    """Per-anchor quality indicator at one test point: mean SNR for MK8000 (exposed by the module);
    for DWM1001C, which exposes no SNR, the consistency of the repeated ranges (minus their std)."""
    q = {}
    for i, a in enumerate(ANCHORS):
        v = g[f"snr_AN{i}"].mean() if system == "MK8000" else -g[f"h_{a}"].std(ddof=1)
        q[a] = -np.inf if np.isnan(v) else v
    return q


# --------------------------------------------------------------------------- pipeline
def locate(system: str, g: pd.DataFrame, estimator: str):
    """Estimated (x, y) and anchors used for the 30 cycles of one test point."""
    h = g[[f"h_{a}" for a in ANCHORS]].mean()              # step 2 (NaN-aware)
    rng = {a: h[f"h_{a}"] for a in ANCHORS if not np.isnan(h[f"h_{a}"])}
    ref, active = select_anchors(rng, quality(system, g))   # step 3
    if estimator == "OLS":
        xy = ols(active, rng, ref)
    elif estimator == "WLS":
        xy = wls(active, rng, ref, {a: g[f"h_{a}"].var(ddof=1) for a in active})
    elif estimator == "RANSAC":
        xy, active = ransac(active, rng, tau(system))
    else:
        raise ValueError(estimator)
    return xy, active


def run(system: str) -> pd.DataFrame:
    """Position error of every test point for the three estimators (long format)."""
    rows = []
    for pid, g in load(system).groupby("point_id"):
        gx, gy, region = g.gt_x.iloc[0], g.gt_y.iloc[0], g.region.iloc[0]
        for est in ESTIMATORS:
            (x, y), used = locate(system, g, est)
            rows.append(dict(system=system, estimator=est, point_id=pid, region=region, gt_x=gx, gt_y=gy,
                             est_x=x, est_y=y,
                             error_cm=100 * np.hypot(x - gx, y - gy), anchors="".join(a[1] for a in used)))
    return pd.DataFrame(rows)


def summarize(errors: pd.DataFrame) -> pd.DataFrame:
    """Error statistics [cm] per system and estimator: mean, sample std and 95th percentile over the 64
    interior points, and the mean over all 68 points."""
    rows = []
    for (system, est), g in errors.groupby(["system", "estimator"], sort=False):
        inner = g[g.region == INTERIOR].error_cm
        rows.append(dict(system=system, estimator=est, mean_64=inner.mean(), std_64=inner.std(ddof=1),
                         p95_64=inner.quantile(0.95), mean_68=g.error_cm.mean()))
    return pd.DataFrame(rows).round(1)


def _error_model(system: str) -> pd.DataFrame:
    """d_measured = a1 * d_real + a0 per anchor, over all test points where the anchor was measured;
    eps = standard deviation of the residuals [m] (unrounded)."""
    df, rows = load(system), []
    for a, (ax, ay, _) in ANCHORS.items():
        g = df.groupby("point_id").agg(h=(f"h_{a}", "mean"), x=("gt_x", "first"), y=("gt_y", "first")).dropna()
        d_real = np.hypot(g.x - ax, g.y - ay)
        a1, a0 = np.polyfit(d_real, g.h, 1)
        rows.append(dict(anchor=a, a0=a0, a1=a1, eps=np.std(g.h - (a1 * d_real + a0), ddof=1), n=len(g)))
    return pd.DataFrame(rows)


def error_model(system: str) -> pd.DataFrame:
    """Per-anchor ranging error model (paper Table II), rounded as printed."""
    return _error_model(system).round(2)


@lru_cache(maxsize=None)
def tau(system: str, k: float = 2.5) -> float:
    """RANSAC consensus threshold [m]: k times the median ranging residual (paper Sec. III-D)."""
    return k * float(np.median(_error_model(system).eps))


# --------------------------------------------------------------------------- command line
def main() -> int:
    errors = pd.concat([run(s) for s in FILES])
    table = summarize(errors)
    RESULTS.mkdir(exist_ok=True)
    errors.round(4).to_csv(RESULTS / "errors_per_point.csv", index=False)
    table.to_csv(RESULTS / "table_iii.csv", index=False)

    paper = [PAPER_TABLE_III[(r.system, r.estimator)] for r in table.itertuples()]
    table["paper_mean_64"], table["paper_p95_64"], table["paper_mean_68"] = zip(*paper)
    table["match"] = ((table.mean_64 == table.paper_mean_64) & (table.p95_64 == table.paper_p95_64)
                      & (table.mean_68 == table.paper_mean_68))
    print("Table III (cm)")
    print(table.to_string(index=False))

    model = pd.concat({s: error_model(s).set_index("anchor") for s in FILES}, names=["system", "anchor"])
    model.to_csv(RESULTS / "table_ii.csv")
    ok2 = []
    for (system, anchor), r in model.iterrows():
        a0, a1, eps, n = PAPER_TABLE_II[(system, anchor)]
        ok2.append(np.isclose([r.a0, r.a1, r.n], [a0, a1, n]).all() and round(abs(r.eps - eps), 2) <= 0.03)
    print(f"\nTable II: a0, a1, n identical and eps within 0.03 m in {sum(ok2)} of {len(ok2)} cells")
    print("RANSAC thresholds [m]: " + ", ".join(f"{s} {tau(s):.3f}" for s in FILES))

    ok = bool(table.match.all()) and all(ok2)
    print("\nMatches the paper (Table III mean / p95 / 68-point mean, Table II): " + ("YES" if ok else "NO"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())

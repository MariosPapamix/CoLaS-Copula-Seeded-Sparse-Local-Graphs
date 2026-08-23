"""
make_demo_data.py -- reproducible semi-synthetic benchmark in exact SNAP
Brightkite file format, generated from the conditional model so that
ground truth is known.  Used by the paper's pipeline-validation section;
the real Brightkite fit is the same pipeline pointed at the SNAP files.

Truth: lambda=2, R=2 km, b=0.25 (psi = 0.0625), same parameters in both
metropolitan regions (so held-out parameter transfer should succeed).

Writes demo_edges.txt.gz, demo_checkins.txt.gz next to this script.
"""
import os, sys, gzip
import numpy as np

BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, BASE)
from brightkite_pipeline import spatial_score, sample_marks_cond
from scipy.spatial import cKDTree

rng = np.random.default_rng(20260812)
TRUE = dict(lam=2.0, Rkm=2.0, b=0.25)

# clustered homes in two metro boxes (SF Bay fit region, Denver held-out)
def metro(lat_c, lon_c, n, cores=3, spread=0.19):
    centers = np.stack([lat_c + rng.normal(0, 0.08, cores),
                        lon_c + rng.normal(0, 0.10, cores)], 1)
    w = rng.dirichlet(np.ones(cores))
    which = rng.choice(cores, n, p=w)
    return centers[which] + rng.normal(0, spread * 0.55, (n, 2))

regions = [("SF", 37.62, -122.17, 3000), ("DEN", 39.80, -105.00, 2200)]
homes, tags = [], []
for tag, la, lo, n in regions:
    homes.append(metro(la, lo, n)); tags += [tag] * n
homes = np.vstack(homes); n_all = homes.shape[0]

edges = []
V_all = np.zeros(n_all)
for tag, la0, la1, lo0, lo1 in (("SF", 37.15, 38.10, -122.75, -121.60),
                                ("DEN", 39.40, 40.20, -105.40, -104.50)):
    m = ((homes[:, 0] >= la0) & (homes[:, 0] <= la1)
         & (homes[:, 1] >= lo0) & (homes[:, 1] <= lo1))
    idx = np.where(m)[0]
    lat_c = np.radians((la0 + la1) / 2)
    X = np.stack([(homes[m, 1] - lo0) * 111.320 * np.cos(lat_c),
                  (homes[m, 0] - la0) * 110.574], 1)
    phi, _ = spatial_score(X)
    V = sample_marks_cond(phi, TRUE["b"], rng)
    V_all[idx] = V
    tree = cKDTree(X)
    P = np.array(sorted(tree.query_pairs(TRUE["Rkm"])))
    p = -np.expm1(-TRUE["lam"] * V[P[:, 0]] * V[P[:, 1]])
    keep = rng.random(p.size) < p
    edges.append(idx[P[keep]])
    print(f"{tag}: n={idx.size} edges={int(keep.sum())} "
          f"mean_deg={2*keep.sum()/idx.size:.2f}")
E = np.vstack(edges)

# activity counts rank-monotone in V, >= 6 so no user is filtered
act = 6 + np.round(V_all * 74).astype(int)
with gzip.open(os.path.join(BASE, "demo_checkins.txt.gz"), "wt") as fh:
    for u in range(n_all):
        for _ in range(act[u]):
            la = homes[u, 0] + rng.normal(0, 0.006)
            lo = homes[u, 1] + rng.normal(0, 0.006)
            fh.write(f"{u}\t2010-01-01T00:00:00Z\t{la:.6f}\t{lo:.6f}\t0\n")
with gzip.open(os.path.join(BASE, "demo_edges.txt.gz"), "wt") as fh:
    for a, b in E:
        fh.write(f"{a}\t{b}\n"); fh.write(f"{b}\t{a}\n")
print("wrote demo_edges.txt.gz, demo_checkins.txt.gz; truth:", TRUE,
      "(psi = %.4f)" % TRUE["b"] ** 2)

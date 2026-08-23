"""Distance-profile null for SF Bay and LA: random pairs matched to the
observed edge-length distribution (length-decile matched); reports C and r.
Appends 'dist_null' entries into brightkite_cities.json."""
import sys, json
import os as _os0
sys.path.insert(0, _os0.path.dirname(_os0.path.abspath(__file__)))
import os as _os
_BASE = _os.environ.get("COLAS_BASE") or _os.path.dirname(
    _os.path.dirname(_os.path.abspath(__file__)))
def _p(*parts):
    q = _os.path.join(_BASE, *parts)
    _os.makedirs(_os.path.dirname(q), exist_ok=True) if parts[:1] in (("runs",), ("figures",), ("revision",)) else None
    return q

import numpy as np
from scipy.spatial import cKDTree
from brightkite_pipeline import load_homes_csv, region_graph, graph_stats

rng = np.random.default_rng(7)
users, lats, lons, counts, E = load_homes_csv(
    _p("data", "loc-brightkite_edges.txt.gz"),
    _p("data", "brightkite_homes.csv.gz"))
CT = json.load(open(_p("runs", "brightkite_cities.json")))
for name in ("SF Bay", "LA"):
    g = region_graph(name, users, lats, lons, counts, E)
    X = g["X"]; n = g["n"]
    d_obs = np.linalg.norm(X[g["R"]] - X[g["C"]], axis=1)
    tree = cKDTree(X)
    outs = []
    for rep in range(12):
        R2, C2 = [], []
        eset = set()
        for ell in d_obs:
            for _ in range(200):
                i = rng.integers(0, n)
                lo, hi = 0.9 * ell, 1.1 * ell + 1e-6
                cand = tree.query_ball_point(X[i], hi)
                cand = [j for j in cand if j != i and np.linalg.norm(X[j]-X[i]) >= lo]
                if cand:
                    j = int(rng.choice(cand))
                    a, b = min(i, j), max(i, j)
                    if (a << 20 | b) not in eset:
                        eset.add(a << 20 | b)
                        R2.append(a); C2.append(b)
                        break
        st = graph_stats(n, np.array(R2), np.array(C2))
        outs.append([st["C"], st["r"]])
        eset.clear()
    outs = np.array(outs)
    CT[name]["dist_null"] = dict(C=float(outs[:,0].mean()), C_sd=float(outs[:,0].std()),
                                 r=float(outs[:,1].mean()), r_sd=float(outs[:,1].std()))
    print(f"{name}: distance-matched null C={outs[:,0].mean():.4f}({outs[:,0].std():.4f}) "
          f"r={outs[:,1].mean():.4f}({outs[:,1].std():.4f})", flush=True)
json.dump(CT, open(_p("runs", "brightkite_cities.json"), "w"), indent=1)
print("updated brightkite_cities.json")

"""Design functionals for Proposition S5.4 (signed-b first-order
information): per city, the mean score elevation of within-range pairs
d1bar = D1(R-hat)/N_pairs with its SE, S_a(lambda-hat), and the
measured Delta_sort(+/-b) asymmetry.  Writes runs/design_functionals.json.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import brightkite_pipeline as bp
from scipy.integrate import dblquad
from scipy.spatial import cKDTree

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

def S_a(lam):
    f = lambda v, u: (1 - np.exp(-lam * u * v)) * np.sqrt(3) * (2 * u - 1)
    val, _ = dblquad(f, 0, 1, 0, 1, epsabs=1e-12)
    return val

def main():
    r8 = {}
    for f in ("cities_r8_brightkite.json", "cities_r8_gowalla.json"):
        r8.update(json.load(open(os.path.join(BASE, "runs", f))))
    out = {}
    plats = [
        ("Brightkite", "data/loc-brightkite_edges.txt.gz",
         "data/brightkite_homes.csv.gz",
         {"SF Bay": (37.15, 38.10, -122.75, -121.60),
          "Denver": (39.40, 40.20, -105.40, -104.50),
          "London": (51.20, 51.75, -0.60, 0.40),
          "LA": (33.60, 34.40, -118.70, -117.60),
          "NYC": (40.45, 41.00, -74.35, -73.55)}),
        ("Gowalla", "data/loc-gowalla_edges.txt.gz",
         "data/gowalla_homes.csv.gz",
         {"Austin": (30.05, 30.60, -98.05, -97.35),
          "Stockholm": (59.20, 59.50, 17.75, 18.35),
          "SF Bay": (37.15, 38.10, -122.75, -121.60),
          "Dallas": (32.55, 33.10, -97.15, -96.45),
          "NYC": (40.45, 41.00, -74.35, -73.55),
          "London": (51.20, 51.75, -0.60, 0.40)})]
    for plat, ef, hf, regions in plats:
        users, lats, lons, counts, E = bp.load_homes_csv(
            os.path.join(BASE, ef), os.path.join(BASE, hf))
        bp.REGIONS.update(regions)
        for name in regions:
            key = f"{plat}/{name}"
            c = r8[key]
            g = bp.region_graph(name, users, lats, lons, counts, E)
            phi, bw = bp.spatial_score(g["X"])
            tree = cKDTree(g["X"])
            Rhat = max(c["theta"][1], 0.1)
            P = np.array(sorted(tree.query_pairs(Rhat)))
            if P.size == 0:
                continue
            s_pair = phi[P[:, 0]] + phi[P[:, 1]]
            out[key] = dict(
                n=int(g["n"]), Rhat=Rhat, n_pairs=int(P.shape[0]),
                d1bar=float(s_pair.mean()),
                d1bar_se=float(s_pair.std() / np.sqrt(s_pair.shape[0])),
                D1_per_node=float(s_pair.sum() / g["n"]),
                S_a_lam=S_a(c["theta"][0]),
                dsort_asym=float(c["dsort_signed"] + c["dsort_flip"]))
            print(key, out[key], flush=True)
    json.dump(out, open(os.path.join(BASE, "runs",
                                     "design_functionals.json"), "w"),
              indent=1)

if __name__ == "__main__":
    main()

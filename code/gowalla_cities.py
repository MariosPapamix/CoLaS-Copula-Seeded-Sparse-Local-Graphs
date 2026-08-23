"""Gowalla multi-city analysis: same per-metro battery as Brightkite
(brightkite_cities.py), on prespecified Gowalla metros including Austin
(Gowalla's home city) and Dallas.  Requires:
  data/loc-gowalla_edges.txt.gz     (uploaded)
  data/gowalla_homes.csv.gz         (from summarize_checkins.py)
Writes runs/gowalla_cities.json
"""
import sys, os, json, time
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
import brightkite_pipeline as bp
from brightkite_cities_lib import run_city   # shared per-city battery

EDGES = _p("data", "loc-gowalla_edges.txt.gz")
HOMES = _p("data", "gowalla_homes.csv.gz")
if not os.path.exists(HOMES):
    sys.exit("gowalla_homes.csv.gz not found - run summarize_checkins.py "
             "on loc-gowalla_totalCheckins.txt locally and upload it.")

bp.REGIONS = {
    "Austin":  (30.05, 30.60, -98.05, -97.35),
    "Stockholm": (59.20, 59.50, 17.75, 18.35),
    "SF Bay":  (37.15, 38.10, -122.75, -121.60),
    "Dallas":  (32.55, 33.10, -97.15, -96.45),
    "NYC":     (40.45, 41.00,  -74.35,  -73.55),
    "London":  (51.20, 51.75,   -0.60,    0.40),
}

t0 = time.time()
users, lats, lons, counts, E = bp.load_homes_csv(EDGES, HOMES)
results = {}
for name in bp.REGIONS:
    out = run_city(name, users, lats, lons, counts, E,
                   rng=np.random.default_rng(20260814), boot=24)
    if out is not None:
        results[name] = out
        print("[%6.0fs] %s done" % (time.time() - t0, name), flush=True)
    with open(_p("runs", "gowalla_cities.json"), "w") as fh:
        json.dump(results, fh, indent=1)
with open(_p("runs", "gowalla_cities.json"), "w") as fh:
    json.dump(results, fh, indent=1)
print("DONE -> gowalla_cities.json")

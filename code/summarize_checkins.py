#!/usr/bin/env python3
"""Pure-stdlib alternative to summarize_checkins.R (no packages needed).

  cd /Users/mariospapamichalis/Downloads
  python3 summarize_checkins.py loc-brightkite_totalCheckins.txt

Accepts the .txt or the original .txt.gz.  Writes brightkite_homes.csv.gz
with columns user, med_lat, med_lon, n_valid, n_checkins.
"""
import sys, gzip, csv
from collections import defaultdict
from statistics import median

import os, re
src = sys.argv[1] if len(sys.argv) > 1 else "loc-brightkite_totalCheckins.txt"
if len(sys.argv) > 2:
    out_name = sys.argv[2]
else:
    base = os.path.basename(src)
    m = re.match(r"loc[-_]?([A-Za-z]+)", base)
    out_name = (m.group(1).lower() if m else "network") + "_homes.csv.gz"
opener = gzip.open if src.endswith(".gz") else open

lats = defaultdict(list); lons = defaultdict(list); n_all = defaultdict(int)
with opener(src, "rt") as fh:
    for line in fh:
        p = line.rstrip("\n").split("\t")
        if len(p) < 5:
            continue
        try:
            u = int(p[0]); la = float(p[2]); lo = float(p[3])
        except ValueError:
            continue
        n_all[u] += 1
        if (la == 0.0 and lo == 0.0) or not (-90 <= la <= 90 and -180 <= lo <= 180):
            continue
        lats[u].append(la); lons[u].append(lo)

with gzip.open(out_name, "wt", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["user", "med_lat", "med_lon", "n_valid", "n_checkins"])
    for u in sorted(lats):
        w.writerow([u, f"{median(lats[u]):.6f}", f"{median(lons[u]):.6f}",
                    len(lats[u]), n_all[u]])
print(f"wrote {out_name}: {len(lats)} users")

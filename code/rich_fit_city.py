"""Fit the richer class to one city, with GOF-before-attribution and
within-class counterfactual attribution.  Usage:
    python3 rich_fit_city.py "SF Bay" brightkite
    python3 rich_fit_city.py "Austin" gowalla
Writes runs/rich_<platform>_<city>.json
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.stats import rankdata
import brightkite_pipeline as bp
from richer_model import RichModel, TARGET_KEYS, fit_rich, gof, rel_resid

CITY = sys.argv[1] if len(sys.argv) > 1 else "SF Bay"
PLAT = sys.argv[2] if len(sys.argv) > 2 else "brightkite"
BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

REGIONS = {
    "brightkite": {"SF Bay": (37.15, 38.10, -122.75, -121.60),
                   "Denver": (39.40, 40.20, -105.40, -104.50),
                   "London": (51.20, 51.75, -0.60, 0.40),
                   "LA":     (33.60, 34.40, -118.70, -117.60),
                   "NYC":    (40.45, 41.00, -74.35, -73.55)},
    "gowalla":    {"Austin":  (30.05, 30.60, -98.05, -97.35),
                   "Stockholm": (59.20, 59.50, 17.75, 18.35),
                   "SF Bay": (37.15, 38.10, -122.75, -121.60),
                   "Dallas": (32.55, 33.10, -97.15, -96.45),
                   "NYC":    (40.45, 41.00, -74.35, -73.55),
                   "London": (51.20, 51.75, -0.60, 0.40)},
}
DATA = {"brightkite": ("data/loc-brightkite_edges.txt.gz",
                       "data/brightkite_homes.csv.gz"),
        "gowalla":    ("data/loc-gowalla_edges.txt.gz",
                       "data/gowalla_homes.csv.gz")}

def main():
    ef, hf = DATA[PLAT]
    users, lats, lons, counts, E = bp.load_homes_csv(
        os.path.join(BASE, ef), os.path.join(BASE, hf))
    bp.REGIONS.update(REGIONS[PLAT])
    g = bp.region_graph(CITY, users, lats, lons, counts, E)
    st = bp.graph_stats(g["n"], g["R"], g["C"])
    phi, bw = bp.spatial_score(g["X"])
    b_dir, se_b, psi_d, se_pd = bp.direct_psi(g["V"], phi)
    d = np.linalg.norm(g["X"][g["R"]] - g["X"][g["C"]], axis=1)
    deg = np.bincount(g["R"], minlength=g["n"]) + np.bincount(g["C"], minlength=g["n"])
    q = np.percentile(d, [25, 50, 75, 90])
    target = [st["ell"], st["C"], st["r"], q[0], q[1], q[2],
              float(deg.std() / deg.mean()), float(np.percentile(deg, 99))]
    rk = rankdata(deg, method="average")
    d1 = np.concatenate([rk[g["R"]], rk[g["C"]]])
    d2 = np.concatenate([rk[g["C"]], rk[g["R"]]])
    obs_spear = float(np.corrcoef(d1, d2)[0, 1])
    obs_held = dict(q90=float(q[3]),
                    hub_share=float(np.sort(deg)[::-1][:max(1, g["n"]//100)].sum()
                                    / max(deg.sum(), 1)),
                    max_deg=int(deg.max()), spearman=obs_spear)
    print(f"[rich] {PLAT}/{CITY} n={g['n']} targets="
          f"{[round(t,3) for t in target]}", flush=True)

    M = RichModel(g["X"], phi)
    t0 = time.time()
    x0 = None
    warm = os.environ.get("RICH_WARM")
    if warm and os.path.exists(warm):
        w0 = json.load(open(warm))["theta"]
        x0 = (w0["lam"], w0["R"], w0["w"], w0["L"], w0["sigma"], w0["tau"], w0["b"])
        print(f"[rich] warm start from {warm}", flush=True)
    fb = os.environ.get("RICH_FIXB")
    th = fit_rich(M, target, sweeps=int(os.environ.get("RICH_SWEEPS", "4")),
                  n_rep=int(os.environ.get("RICH_REPS", "6")), x0=x0,
                  fix_b=(float(fb) if fb not in (None, "", "dir") else
                         (float(b_dir) if fb == "dir" else None)))
    print(f"[rich] fit done in {time.time()-t0:.0f}s", flush=True)

    G = gof(M, th, target, B=100)
    print(f"[rich] GOF max|z| = {G['max_abs_z']:.2f}  z="
          f"{[round(z,1) for z in G['z']]}", flush=True)

    # held-out checks from full sims at theta-hat
    rng = np.random.default_rng(7)
    fulls = M.sim(tuple(th), rng, 30, full=True)
    held = dict(q90=[], hub_share=[], max_deg=[], spearman=[])
    annd_obs, annd_mod = annd_curves(g, deg, fulls, M)
    for f in fulls:
        dg = f["deg"]
        ds = np.sort(dg)[::-1]
        held["hub_share"].append(float(ds[:max(1, M.n//100)].sum() / max(dg.sum(), 1)))
        held["max_deg"].append(float(dg.max()))
        rkm = rankdata(dg, method="average")
        A = f["A"].tocoo()
        upper = A.row < A.col
        rr, cc = A.row[upper], A.col[upper]
        e1 = np.concatenate([rkm[rr], rkm[cc]]); e2 = np.concatenate([rkm[cc], rkm[rr]])
        held["spearman"].append(float(np.corrcoef(e1, e2)[0, 1]) if e1.std() > 0 else 0.0)
    # q90 of edge lengths per sim
    rng2 = np.random.default_rng(8)
    for _ in range(15):
        m_full = M.sim(tuple(th), rng2, 1, full=True)
        if m_full:
            pass
    held_summary = {k: [float(np.mean(v)), float(np.std(v))]
                    for k, v in held.items() if v}

    # within-class counterfactual attribution (each with lam re-matched to ell)
    def recal_lam(theta):
        theta = list(theta)
        for _ in range(8):
            m = M.sim(tuple(theta), np.random.default_rng(999), 6)
            f = target[0] / max(m.mean(0)[0], 1e-9)
            if abs(f - 1) < 0.01:
                break
            theta[0] = float(np.clip(theta[0] * f, 1e-3, 500))
        return theta

    def r_at(theta, reps=40):
        m = M.sim(tuple(theta), np.random.default_rng(31), reps)
        return float(m.mean(0)[2]), float(m.std(0)[2] / np.sqrt(len(m)))

    r_hat, r_hat_se = r_at(th)
    lam, R, w, L, sigma, tau, b = th
    cf = {}
    for name, theta_cf in [
            ("no_closure",   (lam, R, w, L, sigma, 0.0, b)),
            ("no_hubs",      (lam, R, w, L, 0.05, tau, b)),
            ("no_longrange", (lam, R, 0.0, L, sigma, tau, b)),
            ("b_zero",       (lam, R, w, L, sigma, tau, 0.0)),
            ("b_direct",     (lam, R, w, L, sigma, tau, float(b_dir))),
            ("b_direct_flip",(lam, R, w, L, sigma, tau, -float(b_dir)))]:
        tcf = recal_lam(theta_cf)
        r_cf, r_cf_se = r_at(tcf)
        cf[name] = dict(r=r_cf, se=r_cf_se, delta=r_hat - r_cf,
                        theta=[float(x) for x in tcf])
        print(f"[rich] cf {name}: r={r_cf:+.4f} delta={r_hat-r_cf:+.4f}",
              flush=True)

    # graph-information profile in signed b (others at theta-hat, lam recal)
    prof = []
    for bb in np.linspace(-0.3, 0.3, 13):
        tcf = recal_lam((lam, R, w, L, sigma, tau, float(bb)))
        r_b, se_bp = r_at(tcf, reps=25)
        prof.append([float(bb), r_b, se_bp])
    out = dict(city=CITY, platform=PLAT, n=int(g["n"]),
               target=dict(zip(TARGET_KEYS, [float(t) for t in target])),
               obs_held=obs_held,
               theta=dict(zip(["lam", "R", "w", "L", "sigma", "tau", "b"],
                              [float(x) for x in th])),
               b_direct=float(b_dir), se_b_direct=float(se_b),
               gof=G, held=held_summary, r_hat=r_hat, r_hat_se=r_hat_se,
               counterfactuals=cf, b_profile=prof,
               annd=dict(obs=annd_obs, model=annd_mod))
    path = os.path.join(BASE, "runs",
                        f"rich_{PLAT}_{CITY.replace(' ', '')}.json")
    json.dump(out, open(path, "w"), indent=1)
    print(f"[rich] wrote {path}", flush=True)

def annd_curves(g, deg_obs, fulls, M):
    """Mean neighbor degree vs degree, observed and model-averaged."""
    import numpy as np
    def curve(deg, R_idx, C_idx, kmax=16):
        s = np.zeros(max(int(deg.max()), kmax) + 2); c = np.zeros_like(s)
        for a, b in ((R_idx, C_idx), (C_idx, R_idx)):
            for i in range(a.size):
                k = int(deg[a[i]])
                if k <= kmax:
                    s[k] += deg[b[i]]; c[k] += 1
        with np.errstate(invalid="ignore", divide="ignore"):
            return [float(s[k] / c[k]) if c[k] > 0 else None
                    for k in range(1, kmax + 1)]
    obs = curve(deg_obs.astype(float), g["R"], g["C"])
    mods = []
    for f in fulls[:15]:
        A = f["A"].tocoo(); upper = A.row < A.col
        mods.append(curve(f["deg"], A.row[upper], A.col[upper]))
    model = [float(np.mean([m[k] for m in mods if m[k] is not None]))
             if any(m[k] is not None for m in mods) else None
             for k in range(16)]
    return obs, model

if __name__ == "__main__":
    main()

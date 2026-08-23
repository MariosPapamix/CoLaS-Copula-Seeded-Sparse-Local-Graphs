"""City goodness-of-fit battery: identical fit protocol and seeds to
brightkite_cities.py / gowalla_cities.py, with large simulation budgets
(bootstrap draws 200; channels 150, config null 60, Jacobian 60,
fitted/geo moments 120) and goodness-of-fit and counterfactual
outputs per city:
  frac_outside   share of observed edges farther than the fitted radius
  r_resid        observed r minus fitted-model r
  gof_z          r_resid over the bootstrap sd of realized model r at theta-hat
  rails          which parameters sit at optimizer bounds
  dsort_signed   model-implied change in r from b = b_dir (signed), after
                 recalibrating (lambda, R) to observed (ell, C)
  dsort_flip     same with the opposite sign of b (linear-in-b diagnostic)
Writes runs/cities_r8.json incrementally.
"""
import sys, os, json, time
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
import brightkite_pipeline as bp
from brightkite_pipeline import (region_graph, graph_stats, spatial_score,
                                 CondModel, fit_cond, direct_psi)
from brightkite_cities_lib import spearman_r, config_null, model_channels, jac_diag

B_BOOT = int(os.environ.get("R8_BOOT", "150"))
PLATFORM = sys.argv[1] if len(sys.argv) > 1 else None  # "Brightkite"/"Gowalla"/None
_suffix = f"_{PLATFORM.lower()}" if PLATFORM else ""
OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "runs", f"cities_r8{_suffix}.json")

BK_REGIONS = {
    "SF Bay": (37.15, 38.10, -122.75, -121.60),
    "Denver": (39.40, 40.20, -105.40, -104.50),
    "London": (51.20, 51.75, -0.60, 0.40),
    "LA":     (33.60, 34.40, -118.70, -117.60),
    "NYC":    (40.45, 41.00, -74.35, -73.55),
}
GW_REGIONS = {  # cheapest first so partial output accrues early
    "London":  (51.20, 51.75, -0.60, 0.40),
    "NYC":     (40.45, 41.00, -74.35, -73.55),
    "Dallas":  (32.55, 33.10, -97.15, -96.45),
    "SF Bay":  (37.15, 38.10, -122.75, -121.60),
    "Austin":  (30.05, 30.60, -98.05, -97.35),
    "Stockholm": (59.20, 59.50, 17.75, 18.35),
}

def run_city_r8(name, users, lats, lons, counts, E, rng, boot=B_BOOT):
    g = region_graph(name, users, lats, lons, counts, E)
    if g["n"] < 300 or g["R"].size < 200:
        return None
    st = graph_stats(g["n"], g["R"], g["C"])
    deg = np.bincount(g["R"], minlength=g["n"]) + np.bincount(g["C"], minlength=g["n"])
    r_s = spearman_r(deg, g["R"], g["C"])
    (Cn, rn), (Cns, rns) = config_null(g["n"], g["R"], g["C"], rng, n_rep=60)
    phi, bw = spatial_score(g["X"])
    b_h, se_b, psi_d, se_pd = direct_psi(g["V"], phi)
    model = CondModel(g["X"], phi)
    target = np.array([st["ell"], st["C"], st["r"]])
    # identical point-fit protocol (same seeds) as the frozen city runs
    th, q = fit_cond(model, target, np.array([1.0, 1.0, max(b_h, 0.02)]),
                     np.random.default_rng(1))
    th0, _ = fit_cond(model, target, th, np.random.default_rng(2), fix_b=0.0)
    m_fit = model.sim(*th, np.random.default_rng(3), 120).mean(0)
    m_geo = model.sim(*th0, np.random.default_rng(3), 120).mean(0)
    # counterfactual: recalibrate (lam, R) to (ell, C) at the SIGNED direct b
    bd = float(b_h)
    thd, _ = fit_cond(model, target, th0, np.random.default_rng(5), fix_b=bd)
    thf, _ = fit_cond(model, target, th0, np.random.default_rng(5), fix_b=-bd)
    m_d = model.sim(*thd, np.random.default_rng(3), 120).mean(0)
    m_f = model.sim(*thf, np.random.default_rng(3), 120).mean(0)
    # bootstrap of the graph-only fit at theta-hat
    psis, r_draws = [], []
    for _ in range(boot):
        mb = model.sim(*th, rng, 1)[0]
        r_draws.append(float(mb[2]))
        thb, _ = fit_cond(model, mb, th, rng, n_rep=10)
        psis.append(thb[2] ** 2)
    psis = np.array(psis); r_draws = np.array(r_draws)
    ch_mean, ch_sd = model_channels(model, th, rng, n_rep=120)
    jd = jac_diag(model, th, n_rep=60)
    # goodness of fit and rails
    d_edge = np.sqrt(((g["X"][g["R"]] - g["X"][g["C"]]) ** 2).sum(1))
    frac_out = float((d_edge > th[1]).mean())
    r_resid = float(st["r"] - m_fit[2])
    gof_z = float(r_resid / max(r_draws.std(), 1e-9))
    rails = dict(b_low=bool(th[2] <= 1e-9), b_high=bool(th[2] >= 0.55 - 1e-9),
                 R_low=bool(th[1] <= 0.10 + 1e-9),
                 R_high=bool(th[1] >= model.R_max - 1e-9))
    out = dict(n=st["n"], edges=st["edges"], ell=st["ell"], C=st["C"],
               r=st["r"], r_spearman=r_s,
               config_null=dict(C=float(Cn), C_sd=float(Cns),
                                r=float(rn), r_sd=float(rns)),
               b_direct=float(b_h), se_b=float(se_b),
               psi_direct=float(psi_d), se_psi=float(se_pd),
               theta=[float(x) for x in th], psi_graph=float(th[2] ** 2),
               psi_boot_q=[float(np.quantile(psis, x)) for x in (0.025, 0.975)],
               boot_draws=int(psis.size),
               model_moments=[float(x) for x in m_fit],
               r_geo=float(m_geo[2]),
               channels=dict(int_=float(ch_mean[0]), overlap=float(ch_mean[1]),
                             int_sd=float(ch_sd[0]), overlap_sd=float(ch_sd[1])),
               jac=jd, bw_km=float(bw),
               frac_outside=frac_out, med_km=float(np.median(d_edge)),
               r_resid=r_resid, gof_z=gof_z, rails=rails,
               dsort_signed=float(m_d[2] - m_geo[2]),
               dsort_flip=float(m_f[2] - m_geo[2]),
               b_dir_signed=bd)
    print(f"[r8] {name}: boot={psis.size} r_resid={r_resid:+.3f} "
          f"gof_z={gof_z:+.1f} frac_out={frac_out:.3f} "
          f"dsort={out['dsort_signed']:+.4f}/{out['dsort_flip']:+.4f}", flush=True)
    return out

def main():
    results = {}
    if os.path.exists(OUT):
        results = json.load(open(OUT))
    jobs = []
    platforms = [
        ("Brightkite", "data/loc-brightkite_edges.txt.gz",
         "data/brightkite_homes.csv.gz", BK_REGIONS),
        ("Gowalla", "data/loc-gowalla_edges.txt.gz",
         "data/gowalla_homes.csv.gz", GW_REGIONS)]
    if PLATFORM:
        platforms = [p for p in platforms if p[0] == PLATFORM]
    for plat, edges_f, homes_f, regions in platforms:
        base = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        users, lats, lons, counts, E = bp.load_homes_csv(
            os.path.join(base, edges_f), os.path.join(base, homes_f))
        bp.REGIONS.update(regions)
        for name in regions:
            key = f"{plat}/{name}"
            if key in results:
                print(f"[r8] skip {key} (done)", flush=True)
                continue
            t0 = time.time()
            out = run_city_r8(name, users, lats, lons, counts, E,
                              rng=np.random.default_rng(20260815))
            if out is not None:
                results[key] = out
                with open(OUT, "w") as fh:
                    json.dump(results, fh, indent=1)
            print(f"[r8] {key} finished in {time.time()-t0:.0f}s", flush=True)
    print("[r8] ALL DONE", flush=True)

if __name__ == "__main__":
    main()

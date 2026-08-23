"""OpenAlex Netherlands collaboration network: descriptive battery,
direct estimate, and base-class conditional fit (two-view analysis).

Data (data/oa_homes.csv.gz, data/oa_edges.txt.gz) are produced by the
OpenAlex fetch script: journal articles 2019--2023 with a Netherlands
institution in the authorship lineage, each author placed at their
institution coordinates, per-author work counts within the corpus, and
co-authorship edges (works with more than 25 authors skipped).

Positions are institution point masses.  A deterministic jitter
(sigma = JIT_KM km, seed JIT_SEED) breaks the atoms so that
nearest-neighbor bandwidths and pair distances are well defined; the
density bandwidth is floored at BW_FLOOR km.  Atom identity (clustered
standard errors, top-atom share) uses raw coordinates rounded to three
decimals (about 100 m).

Stages
------
  descriptives   per region: graph moments, rank assortativity,
                 configuration null (small regions), edge-length
                 quantiles, direct estimate of b with naive and
                 atom-clustered standard errors, edge density score
                 sum D1-bar (signed-b diagnostic).
  OA_BASE=1      base-class conditional fit (hard disc, lambda V_i V_j)
                 on the seeded Netherlands node subsample, with the
                 goodness-of-fit outputs (frac_outside, r_resid,
                 rails, signed counterfactual Delta-sort).

Usage: python3 oa_analysis.py            (descriptives)
       OA_BASE=1 python3 oa_analysis.py  (descriptives if missing + fit)
Writes runs/oa_descriptives.json, runs/oa_base_NLsub.json
"""
import sys, os, gzip, json, time, csv
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.stats import rankdata
import brightkite_pipeline as bp
from brightkite_cities_lib import spearman_r, config_null

_BASE = os.environ.get("COLAS_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))
RUNS = os.path.join(_BASE, "runs")
os.makedirs(RUNS, exist_ok=True)

JIT_SEED, JIT_KM = 20260823, 0.30
BW_FLOOR = 0.5
SUB_SEED, SUB_M = 101, 6500
SQRT3 = np.sqrt(3.0)

OA_REGIONS = {
    "Netherlands":  (50.70, 53.60, 3.30, 7.30),
    "Randstad":     (51.80, 52.45, 4.25, 5.30),
    "Amsterdam":    (52.28, 52.42, 4.75, 5.00),
    "SouthHolland": (51.85, 52.20, 4.25, 4.75),
    "Utrecht":      (52.02, 52.15, 5.05, 5.30),
    "Groningen":    (53.15, 53.30, 6.50, 6.70),
    "Eindhoven":    (51.38, 51.52, 5.35, 5.65),
}

T0 = time.time()
def log(m): print("[%7.1fs] %s" % (time.time() - T0, m), flush=True)

# ---------- load -------------------------------------------------------
def load_oa():
    hp = os.path.join(_BASE, "data", "oa_homes.csv.gz")
    ep = os.path.join(_BASE, "data", "oa_edges.txt.gz")
    uid, lat0, lon0, cnt = [], [], [], []
    with gzip.open(hp, "rt") as fh:
        for row in csv.DictReader(fh):
            uid.append(int(row["user_id"]))
            lat0.append(float(row["lat"]))
            lon0.append(float(row["lon"]))
            cnt.append(int(row["checkins"]))
    uid = np.array(uid); lat0 = np.array(lat0)
    lon0 = np.array(lon0); cnt = np.array(cnt)
    rng = np.random.default_rng(JIT_SEED)
    lat = lat0 + rng.normal(0.0, JIT_KM, uid.size) / 110.574
    lon = lon0 + rng.normal(0.0, JIT_KM, uid.size) / (
        111.320 * np.cos(np.radians(np.clip(lat0, -80, 80))))
    E = []
    with gzip.open(ep, "rt") as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 2:
                a, b = int(p[0]), int(p[1])
                if a != b:
                    E.append((min(a, b), max(a, b)))
    E = np.unique(np.array(E, dtype=np.int64), axis=0)
    log("authors %d edges %d" % (uid.size, E.shape[0]))
    return uid, lat0, lon0, lat, lon, cnt, E

def region(name, uid, lat0, lon0, lat, lon, cnt, E):
    la0, la1, lo0, lo1 = OA_REGIONS[name]
    m = (lat0 >= la0) & (lat0 <= la1) & (lon0 >= lo0) & (lon0 <= lo1)
    u = uid[m]
    idx = {a: i for i, a in enumerate(u)}
    n = u.size
    lat_c = np.radians((la0 + la1) / 2)
    X = np.stack([(lon[m] - lo0) * 111.320 * np.cos(lat_c),
                  (lat[m] - la0) * 110.574], axis=1)
    atom = np.array(["%.3f_%.3f" % (a, b)
                     for a, b in zip(lat0[m], lon0[m])])
    keep = np.array([(a in idx) and (b in idx) for a, b in E])
    Er = E[keep]
    R = np.array([idx[a] for a, _ in Er], dtype=np.int64)
    C = np.array([idx[b] for _, b in Er], dtype=np.int64)
    V = (rankdata(cnt[m], method="average") - 0.5) / n
    return dict(name=name, n=n, X=X, V=V, R=R, C=C,
                atom=atom, counts=cnt[m])

def subsample(g, m_sub, seed):
    rng = np.random.default_rng(seed)
    sel = np.sort(rng.choice(g["n"], size=m_sub, replace=False))
    pos = -np.ones(g["n"], dtype=np.int64)
    pos[sel] = np.arange(m_sub)
    ok = (pos[g["R"]] >= 0) & (pos[g["C"]] >= 0)
    return dict(name=g["name"] + "-sub", n=m_sub, X=g["X"][sel],
                V=(rankdata(g["counts"][sel], method="average") - 0.5) / m_sub,
                R=pos[g["R"][ok]], C=pos[g["C"][ok]],
                atom=g["atom"][sel], counts=g["counts"][sel])

# ---------- battery ----------------------------------------------------
def phi_floor(X):
    phi, bw = bp.spatial_score(X)
    if bw < BW_FLOOR:
        phi, bw = bp.spatial_score(X, bw=BW_FLOOR)
    return phi, bw

def battery(g, cfg=True):
    st = bp.graph_stats(g["n"], g["R"], g["C"])
    deg = (np.bincount(g["R"], minlength=g["n"])
           + np.bincount(g["C"], minlength=g["n"]))
    phi, bw = phi_floor(g["X"])
    b_h, se_b, psi_d, se_pd = bp.direct_psi(g["V"], phi)
    # atom-clustered standard error for b_hat (nodes at one institution
    # share phi and have dependent counts)
    y = SQRT3 * (2 * g["V"] - 1) * phi
    _, inv = np.unique(g["atom"], return_inverse=True)
    sums = np.bincount(inv, weights=y - y.mean())
    se_cl = float(np.sqrt((sums ** 2).sum()) / g["n"])
    atom_sizes = np.bincount(inv)
    out = dict(n=int(g["n"]), edges=int(st["edges"]), ell=st["ell"],
               C=st["C"], r=st["r"], deg_max=int(st["deg_max"]),
               deg_cv=float(deg.std() / max(deg.mean(), 1e-9)),
               p99=float(np.percentile(deg, 99)),
               n_atoms=int(atom_sizes.size),
               top_atom_share=float(atom_sizes.max() / g["n"]),
               bw_km=float(bw), b_direct=float(b_h), se_b=float(se_b),
               se_b_cluster=se_cl, psi_direct=float(psi_d),
               se_psi=float(se_pd),
               works_median=float(np.median(g["counts"])),
               works_max=int(g["counts"].max()))
    if st["edges"] >= 200:
        out["r_spearman"] = spearman_r(deg, g["R"], g["C"])
        d_edge = np.sqrt(((g["X"][g["R"]] - g["X"][g["C"]]) ** 2).sum(1))
        q = np.percentile(d_edge, [25, 50, 75, 90])
        out["edge_km_q"] = [float(x) for x in q]
        out["frac_within_1km"] = float((d_edge <= 1.0).mean())
        out["frac_beyond_6km"] = float((d_edge > 6.0).mean())
        # signed-b diagnostic D1-bar: mean of phi_i + phi_j over edges
        d1 = phi[g["R"]] + phi[g["C"]]
        out["d1_bar"] = float(d1.mean())
        out["d1_se"] = float(d1.std(ddof=1) / np.sqrt(d1.size))
        if cfg and st["edges"] <= 60000:
            (Cn, rn), (Cns, rns) = config_null(
                g["n"], g["R"], g["C"], np.random.default_rng(11), n_rep=20)
            out["config_null"] = dict(C=float(Cn), C_sd=float(Cns),
                                      r=float(rn), r_sd=float(rns))
    return out

# ---------- base-class conditional fit ---------------------------------
def base_fit(g, out_path):
    st = bp.graph_stats(g["n"], g["R"], g["C"])
    phi, bw = phi_floor(g["X"])
    b_h, se_b, psi_d, se_pd = bp.direct_psi(g["V"], phi)
    target = np.array([st["ell"], st["C"], st["r"]])
    log("base fit: building CondModel (n=%d) ..." % g["n"])
    model = bp.CondModel(g["X"], phi)
    log("candidate pairs %d, R_max %.3f km, triangles %d"
        % (model.P.shape[0], model.R_max, model.T1.size))
    th, q = bp.fit_cond(model, target,
                        np.array([1.0, 1.0, max(b_h, 0.02)]),
                        np.random.default_rng(1))
    th0, _ = bp.fit_cond(model, target, th, np.random.default_rng(2),
                         fix_b=0.0)
    bd = float(b_h)
    thd, _ = bp.fit_cond(model, target, th0, np.random.default_rng(5),
                         fix_b=min(abs(bd), 0.55) * np.sign(bd))
    thf, _ = bp.fit_cond(model, target, th0, np.random.default_rng(5),
                         fix_b=-min(abs(bd), 0.55) * np.sign(bd))
    m_fit = model.sim(*th, np.random.default_rng(3), 120).mean(0)
    m_geo = model.sim(*th0, np.random.default_rng(3), 120).mean(0)
    m_d = model.sim(*thd, np.random.default_rng(3), 120).mean(0)
    m_f = model.sim(*thf, np.random.default_rng(3), 120).mean(0)
    r_draws = np.array([model.sim(*th, np.random.default_rng(100 + i), 1)[0][2]
                        for i in range(60)])
    d_edge = np.sqrt(((g["X"][g["R"]] - g["X"][g["C"]]) ** 2).sum(1))
    out = dict(region=g["name"], n=int(g["n"]), edges=int(st["edges"]),
               target=[float(x) for x in target],
               theta=[float(x) for x in th], psi_graph=float(th[2] ** 2),
               model_moments=[float(x) for x in m_fit],
               r_geo=float(m_geo[2]),
               r_resid=float(st["r"] - m_fit[2]),
               gof_z=float((st["r"] - m_fit[2]) / max(r_draws.std(), 1e-9)),
               frac_outside=float((d_edge > th[1]).mean()),
               med_km=float(np.median(d_edge)),
               rails=dict(b_low=bool(th[2] <= 1e-9),
                          b_high=bool(th[2] >= 0.55 - 1e-9),
                          R_low=bool(th[1] <= 0.10 + 1e-9),
                          R_high=bool(th[1] >= model.R_max - 1e-9)),
               dsort_signed=float(m_d[2] - m_geo[2]),
               dsort_flip=float(m_f[2] - m_geo[2]),
               b_dir_signed=bd, R_max_cand=float(model.R_max))
    json.dump(out, open(out_path, "w"), indent=1)
    log("base fit: theta=%s psi_gr=%.5f r_resid=%+.4f frac_out=%.3f"
        % (np.round(th, 4).tolist(), th[2] ** 2, out["r_resid"],
           out["frac_outside"]))
    return out

# ---------- preprocessing sensitivity ----------------------------------
def sensitivity(loaded):
    """Direct estimate at the country scale across jitter scales,
    bandwidth floors, and a replacement jitter seed (OA_SENS=1)."""
    uid, lat0, lon0, _, _, cnt, E = loaded
    la0, la1, lo0, lo1 = OA_REGIONS["Netherlands"]
    m = (lat0 >= la0) & (lat0 <= la1) & (lon0 >= lo0) & (lon0 <= lo1)
    n = int(m.sum())
    lat_c = np.radians((la0 + la1) / 2)
    V = (rankdata(cnt[m], method="average") - 0.5) / n
    atom = np.array(["%.3f_%.3f" % (a, b) for a, b in
                     zip(lat0[m], lon0[m])])
    _, inv = np.unique(atom, return_inverse=True)
    cfgs = ([(s, 0.5, JIT_SEED) for s in (0.15, 0.30, 0.60)]
            + [(0.30, 1.0, JIT_SEED), (0.30, 2.0, JIT_SEED),
               (0.60, 1.0, JIT_SEED), (0.30, 0.5, 20260901)])
    out = []
    for sj, fl, seed in cfgs:
        rng = np.random.default_rng(seed)
        la = lat0 + rng.normal(0.0, sj, uid.size) / 110.574
        lo = lon0 + rng.normal(0.0, sj, uid.size) / (
            111.320 * np.cos(np.radians(np.clip(lat0, -80, 80))))
        X = np.stack([(lo[m] - lo0) * 111.320 * np.cos(lat_c),
                      (la[m] - la0) * 110.574], axis=1)
        phi, bw = bp.spatial_score(X)
        if bw < fl:
            phi, bw = bp.spatial_score(X, bw=fl)
        b, se, psi, sep = bp.direct_psi(V, phi)
        y = SQRT3 * (2 * V - 1) * phi
        sums = np.bincount(inv, weights=y - y.mean())
        se_cl = float(np.sqrt((sums ** 2).sum()) / n)
        out.append(dict(jit_km=sj, floor_km=fl, seed=int(seed),
                        bw_km=float(bw), b_direct=float(b),
                        se_b_cluster=se_cl))
        log("sens jit=%.2f floor=%.1f seed=%d: b=%+.4f (cl %.4f) "
            "bw=%.2f" % (sj, fl, seed, b, se_cl, bw))
    json.dump(out, open(os.path.join(RUNS, "oa_sensitivity.json"), "w"),
              indent=1)
    log("wrote oa_sensitivity.json")

# ---------- main -------------------------------------------------------
def build_nlsub(loaded=None):
    if loaded is None:
        loaded = load_oa()
    gNL = region("Netherlands", *loaded)
    return subsample(gNL, min(SUB_M, gNL["n"]), SUB_SEED), loaded

def main():
    desc_path = os.path.join(RUNS, "oa_descriptives.json")
    do_base = os.environ.get("OA_BASE") == "1"
    loaded = load_oa()
    if os.environ.get("OA_SENS") == "1":
        sensitivity(loaded)
        return
    if not (do_base and os.path.exists(desc_path)):
        results = dict(meta=dict(jitter_km=JIT_KM, jitter_seed=JIT_SEED,
                                 bw_floor_km=BW_FLOOR, sub_m=SUB_M,
                                 sub_seed=SUB_SEED), regions={})
        for name in OA_REGIONS:
            g = region(name, *loaded)
            results["regions"][name] = battery(g)
            b = results["regions"][name]
            log("%-13s n=%d E=%d ell=%.2f C=%.3f r=%+.3f | b_dir=%+.4f "
                "(se %.4f, cl %.4f) atoms=%d top=%.2f"
                % (name, b["n"], b["edges"], b["ell"], b["C"], b["r"],
                   b["b_direct"], b["se_b"], b["se_b_cluster"],
                   b["n_atoms"], b["top_atom_share"]))
            json.dump(results, open(desc_path, "w"), indent=1)
        gsub, _ = build_nlsub(loaded)
        results["regions"]["NL-sub"] = battery(gsub)
        b = results["regions"]["NL-sub"]
        log("%-13s n=%d E=%d ell=%.2f C=%.3f r=%+.3f | b_dir=%+.4f "
            "(se %.4f, cl %.4f)"
            % ("NL-sub", b["n"], b["edges"], b["ell"], b["C"], b["r"],
               b["b_direct"], b["se_b"], b["se_b_cluster"]))
        json.dump(results, open(desc_path, "w"), indent=1)
        log("wrote " + desc_path)
    if do_base:
        gsub, _ = build_nlsub(loaded)
        base_fit(gsub, os.path.join(RUNS, "oa_base_NLsub.json"))

if __name__ == "__main__":
    main()

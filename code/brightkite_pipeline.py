"""
brightkite_pipeline.py
======================
Empirical application of the conditional framework to the SNAP Brightkite
location-based social network, conditional on observed home locations.

Stages
------
 1. Parse check-ins -> per-user home (median lat/lon), activity count.
 2. Select two dense metro regions (fit region + held-out region).
 3. Build region graphs (homes in box, >=5 check-ins, induced edges).
 4. Direct estimate of dependence psi from (X_i, V_i) pairs alone
    (V = activity rank; spatial score phi = standardized log home-density).
 5. Graph-only fit of (lambda, R_km, b) by simulated moment matching of
    (mean degree, transitivity, Newman r), positions held fixed.
 6. psi=0 ablation with (lambda, R) recalibrated to (ell, C):
    r_geo = mixing produced by geometry + density alone.
    Headline: 100 * r_geo / r_obs  "% of observed assortativity is
    geometric baseline".
 7. Held-out validation: non-fitted statistics on the fit region and
    parameter transfer to the held-out region.

Usage: python3 brightkite_pipeline.py <edges.txt.gz> <checkins.txt.gz>
Writes runs/brightkite_results.json (+ npz arrays).
"""
import sys, gzip, json, time
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
import scipy.sparse as sp

OUT = _p("runs") + _os.sep
RNG = np.random.default_rng(20260812)
SQRT3 = np.sqrt(3.0)

REGIONS = {
    "SF Bay":  (37.15, 38.10, -122.75, -121.60),
    "Denver":  (39.40, 40.20, -105.40, -104.50),
    "London":  (51.20, 51.75,   -0.60,    0.40),
    "LA":      (33.60, 34.40, -118.70, -117.60),
    "NYC":     (40.45, 41.00,  -74.35,  -73.55),
}

def log(m): print("[%7.1fs] %s" % (time.time() - T0, m), flush=True)
T0 = time.time()

# ---------- 1. parse ---------------------------------------------------
def load_homes_csv(edges_path, homes_path, min_checkins=5):
    """Load the per-user homes summary produced by summarize_checkins.R/.py
    (columns user, med_lat, med_lon, n_valid, n_checkins) plus the edges."""
    import csv
    log("loading homes summary ...")
    opener = gzip.open if homes_path.endswith(".gz") else open
    users, lats, lons, counts = [], [], [], []
    with opener(homes_path, "rt") as fh:
        rd = csv.DictReader(fh)
        for row in rd:
            if int(row["n_valid"]) < min_checkins:
                continue
            users.append(int(row["user"]))
            lats.append(float(row["med_lat"]))
            lons.append(float(row["med_lon"]))
            counts.append(int(row["n_checkins"]))
    users = np.array(users); lats = np.array(lats)
    lons = np.array(lons); counts = np.array(counts)
    log("users with home: %d" % users.size)
    log("parsing edges ...")
    E = []
    with gzip.open(edges_path, "rt") as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 2:
                a, b = int(p[0]), int(p[1])
                if a < b:
                    E.append((a, b))
    E = np.array(E, dtype=np.int64)
    log("undirected edges: %d" % E.shape[0])
    return users, lats, lons, counts, E

def parse(edges_path, checkins_path, min_checkins=5):
    log("parsing check-ins ...")
    users, lats, lons = [], [], []
    with gzip.open(checkins_path, "rt") as fh:
        cur, la, lo = None, [], []
        for line in fh:
            p = line.rstrip("\n").split("\t")
            if len(p) < 5:
                continue
            u = int(p[0])
            try:
                a, b = float(p[2]), float(p[3])
            except ValueError:
                continue
            if a == 0.0 and b == 0.0:
                continue
            if not (-90 <= a <= 90 and -180 <= b <= 180):
                continue
            if u != cur:
                if cur is not None and len(la) >= min_checkins:
                    users.append(cur); lats.append(np.median(la)); lons.append(np.median(lo))
                cur, la, lo = u, [], []
            la.append(a); lo.append(b)
        if cur is not None and len(la) >= min_checkins:
            users.append(cur); lats.append(np.median(la)); lons.append(np.median(lo))
    # counts need a second cheap pass (streaming above kept only coords)
    cnt = {}
    with gzip.open(checkins_path, "rt") as fh:
        for line in fh:
            u = line.split("\t", 1)[0]
            try:
                u = int(u)
            except ValueError:
                continue
            cnt[u] = cnt.get(u, 0) + 1
    users = np.array(users); lats = np.array(lats); lons = np.array(lons)
    counts = np.array([cnt.get(u, 0) for u in users])
    log("users with home: %d" % users.size)
    log("parsing edges ...")
    E = []
    with gzip.open(edges_path, "rt") as fh:
        for line in fh:
            p = line.split()
            if len(p) >= 2:
                a, b = int(p[0]), int(p[1])
                if a < b:
                    E.append((a, b))
    E = np.array(E, dtype=np.int64)
    log("undirected edges: %d" % E.shape[0])
    return users, lats, lons, counts, E

# ---------- 2/3. regions ----------------------------------------------
def region_graph(name, users, lats, lons, counts, E):
    la0, la1, lo0, lo1 = REGIONS[name]
    m = (lats >= la0) & (lats <= la1) & (lons >= lo0) & (lons <= lo1)
    uid = users[m]
    idx = {u: i for i, u in enumerate(uid)}
    n = uid.size
    # km coordinates (equirectangular about region center)
    lat_c = np.radians((la0 + la1) / 2)
    X = np.stack([(lons[m] - lo0) * 111.320 * np.cos(lat_c),
                  (lats[m] - la0) * 110.574], axis=1)
    keep = np.array([(a in idx) and (b in idx) for a, b in E])
    Er = E[keep]
    R = np.array([idx[a] for a, b in Er]); C = np.array([idx[b] for a, b in Er])
    from scipy.stats import rankdata
    V = (rankdata(counts[m], method="average") - 0.5) / n  # midrank in (0,1)
    return dict(name=name, n=n, X=X, V=V, R=R, C=C)

# ---------- graph statistics ------------------------------------------
def graph_stats(n, R, C):
    deg = np.bincount(R, minlength=n) + np.bincount(C, minlength=n)
    E_n = R.size
    S2 = 0.5 * float((deg.astype(float) * (deg - 1)).sum())
    A = sp.csr_matrix((np.ones(2 * R.size, dtype=np.float32),
                       (np.concatenate([R, C]), np.concatenate([C, R]))), shape=(n, n))
    T = float((A @ A).multiply(A).sum()) / 6.0
    H3 = float((deg.astype(float) ** 3).sum())
    H11 = float(deg.astype(float) @ (A @ deg.astype(float)))
    d1 = np.concatenate([deg[R], deg[C]]).astype(float)
    d2 = np.concatenate([deg[C], deg[R]]).astype(float)
    r = float(np.corrcoef(d1, d2)[0, 1]) if d1.size > 2 and d1.std() > 0 else 0.0
    C_hat = 3.0 * T / S2 if S2 > 0 else 0.0
    return dict(n=n, edges=E_n, ell=2.0 * E_n / n, C=C_hat, r=r,
                S2=S2, T=T, H3=H3, H11=H11,
                deg_var=float(deg.var()), deg_max=int(deg.max()))

# ---------- 4. spatial score + direct psi ------------------------------
def spatial_score(X, bw=None):
    """standardized log kernel density of homes (leave-one-out-ish)."""
    tree = cKDTree(X)
    if bw is None:
        # median distance to 15th neighbor
        d, _ = tree.query(X, k=16)
        bw = np.median(d[:, -1])
    counts = tree.query_ball_point(X, r=bw, return_length=True) - 1
    dens = np.log1p(counts)
    phi = (dens - dens.mean()) / (dens.std() + 1e-12)
    return phi, bw

def direct_psi(V, phi):
    a = SQRT3 * (2 * V - 1)
    s = a * phi
    b_hat = s.mean()
    se_b = s.std(ddof=1) / np.sqrt(s.size)
    psi_hat = b_hat ** 2
    se_psi = 2 * abs(b_hat) * se_b            # delta method
    return b_hat, se_b, psi_hat, se_psi

# ---------- 5. conditional simulator + fit -----------------------------
def sample_marks_cond(phi, b, rng):
    u = rng.random(phi.size)
    s = SQRT3 * b * phi
    s = np.clip(s, -0.999, 0.999)             # keep density nonnegative
    small = np.abs(s) < 1e-9
    ss = np.where(small, 1.0, s)
    disc = (1 - ss) ** 2 + 4 * ss * u
    v = (-(1 - ss) + np.sqrt(np.maximum(disc, 0))) / (2 * ss)
    return np.where(small, u, np.clip(v, 0, 1))

class CondModel:
    """Conditional-on-positions model, d=2: hard disc radius R_km,
    p = 1 - exp(-lam * V_i V_j) inside the disc.  Fast core: candidate
    pairs and candidate-triangle pair-index triples are precomputed once;
    each simulated graph is then a few boolean vector operations."""
    def __init__(self, X, phi, R_max=6.0, max_cand_deg=60.0):
        self.X = X; self.phi = phi; self.n = X.shape[0]
        self.tree = cKDTree(X)
        P = np.array(sorted(self.tree.query_pairs(R_max)))
        d = np.linalg.norm(X[P[:, 0]] - X[P[:, 1]], axis=1)
        o = np.argsort(d, kind="stable")
        P = P[o]; d = d[o]
        # adaptive cap: keep candidate mean degree bounded
        cap = int(max_cand_deg * self.n / 2)
        if P.shape[0] > cap:
            P = P[:cap]; d = d[:cap]
            R_max = float(d[-1])
        self.R_max = R_max
        self.P = P; self.pd = d
        m = self.P.shape[0]
        # pair-key -> pair index
        key = self.P[:, 0].astype(np.int64) * self.n + self.P[:, 1]
        self.keymap = dict(zip(key.tolist(), range(m)))
        # adjacency lists of the candidate graph
        adj = [[] for _ in range(self.n)]
        for idx in range(m):
            a, b = int(self.P[idx, 0]), int(self.P[idx, 1])
            adj[a].append(b); adj[b].append(a)
        adj = [np.array(sorted(v), dtype=np.int64) for v in adj]
        # enumerate candidate triangles i<j<k as pair-index triples
        t1, t2, t3 = [], [], []
        keyset = set(key.tolist())
        for i in range(self.n):
            nb = adj[i]; nb = nb[nb > i]
            if nb.size < 2:
                continue
            for aa in range(nb.size - 1):
                j = int(nb[aa])
                rest = nb[aa + 1:]
                for k_ in rest:
                    if int(j) * self.n + int(k_) in keyset:
                        t1.append(self.keymap[i * self.n + j])
                        t2.append(self.keymap[i * self.n + int(k_)])
                        t3.append(self.keymap[int(j) * self.n + int(k_)])
        self.T1 = np.array(t1, dtype=np.int64)
        self.T2 = np.array(t2, dtype=np.int64)
        self.T3 = np.array(t3, dtype=np.int64)

    def sim(self, lam, Rkm, b, rng, n_rep=1):
        within = self.pd <= Rkm
        P0, P1 = self.P[:, 0], self.P[:, 1]
        out = []
        for _ in range(n_rep):
            V = sample_marks_cond(self.phi, b, rng)
            p = -np.expm1(-lam * V[P0] * V[P1])
            keep = within & (rng.random(p.size) < p)
            deg = (np.bincount(P0[keep], minlength=self.n)
                   + np.bincount(P1[keep], minlength=self.n)).astype(np.float64)
            E_n = int(keep.sum())
            S2 = 0.5 * float((deg * (deg - 1)).sum())
            T = float(np.count_nonzero(keep[self.T1] & keep[self.T2] & keep[self.T3]))
            d1 = np.concatenate([deg[P0[keep]], deg[P1[keep]]])
            d2 = np.concatenate([deg[P1[keep]], deg[P0[keep]]])
            r = float(np.corrcoef(d1, d2)[0, 1]) if d1.size > 2 and d1.std() > 0 else 0.0
            Ch = 3.0 * T / S2 if S2 > 0 else 0.0
            out.append([2.0 * E_n / self.n, Ch, r])
        return np.array(out)

def fit_cond(model, target, x0, rng, n_rep=12, fix_b=None, tol=2e-3):
    """Nested fit: inner 2-D secant on (lam, R) matching (ell, C) for a
    given b; outer 1-D search on b matching r.  CRN via fixed seed."""
    def inner(b, lam, Rkm):
        """Match (ell, C).  Under independent thinning of the candidate
        geometric graph, C scales with the retention level (lambda) while
        ell scales with retention x candidate degree (lambda and R):
        so lambda tracks C and R tracks ell at the new lambda."""
        for it in range(20):
            m = model.sim(lam, Rkm, b, np.random.default_rng(999), n_rep).mean(0)
            fl = target[0] / max(m[0], 1e-9)
            fc = target[1] / max(m[1], 1e-9)
            if abs(fl - 1) < tol and abs(fc - 1) < tol:
                break
            lam = float(np.clip(lam * fc ** 1.2, 0.02, 80.0))
            m2 = model.sim(lam, Rkm, b, np.random.default_rng(999), n_rep).mean(0)
            fl2 = target[0] / max(m2[0], 1e-9)
            Rkm = float(np.clip(Rkm * fl2 ** 0.55, 0.10, model.R_max))
        m = model.sim(lam, Rkm, b, np.random.default_rng(999), n_rep).mean(0)
        return lam, Rkm, m
    if fix_b is not None:
        lam, Rkm, m = inner(fix_b, x0[0], x0[1])
        d = m - target
        return np.array([lam, Rkm, fix_b]), float(d @ np.diag([1, 400, 400]) @ d)
    # outer: secant on b for r
    b_lo, b_hi = 0.0, 0.55
    lam, Rkm = x0[0], x0[1]
    b = float(np.clip(x0[2], 0.01, 0.5))
    hist = []
    for it in range(14):
        lam, Rkm, m = inner(b, lam, Rkm)
        hist.append((b, m[2]))
        if abs(m[2] - target[2]) < tol:
            break
        if len(hist) >= 2 and abs(hist[-1][1] - hist[-2][1]) > 1e-9:
            (bA, rA), (bB, rB) = hist[-2], hist[-1]
            b_new = bB + (target[2] - rB) * (bB - bA) / (rB - rA)
        else:
            b_new = b + (0.08 if m[2] < target[2] else -0.08)
        b = float(np.clip(b_new, b_lo, b_hi))
    lam, Rkm, m = inner(b, lam, Rkm)
    d = m - target
    return np.array([lam, Rkm, b]), float(d @ np.diag([1, 400, 400]) @ d)

# ---------- main -------------------------------------------------------
def main(edges_path, checkins_path, out_prefix="brightkite"):
    if ".csv" in checkins_path:
        users, lats, lons, counts, E = load_homes_csv(edges_path, checkins_path)
    else:
        users, lats, lons, counts, E = parse(edges_path, checkins_path)
    sizes = {}
    for name in REGIONS:
        la0, la1, lo0, lo1 = REGIONS[name]
        m = (lats >= la0) & (lats <= la1) & (lons >= lo0) & (lons <= lo1)
        sizes[name] = int(m.sum())
    log("region sizes: %s" % sizes)
    ordered = sorted(sizes, key=sizes.get, reverse=True)
    fit_name, held_name = ordered[0], ordered[1]
    results = dict(region_sizes=sizes, fit_region=fit_name, held_region=held_name)

    out_np = {}
    for name in (fit_name, held_name):
        g = region_graph(name, users, lats, lons, counts, E)
        st = graph_stats(g["n"], g["R"], g["C"])
        phi, bw = spatial_score(g["X"])
        b_hat, se_b, psi_d, se_psi = direct_psi(g["V"], phi)
        log(f"{name}: n={g['n']} edges={st['edges']} ell={st['ell']:.3f} "
            f"C={st['C']:.4f} r={st['r']:.4f} | direct b={b_hat:.4f}({se_b:.4f}) "
            f"psi={psi_d:.5f}({se_psi:.5f}) bw={bw:.2f}km")
        results[name] = dict(stats=st, b_direct=b_hat, se_b=se_b,
                             psi_direct=psi_d, se_psi=se_psi, bw_km=bw)
        out_np[f"X_{name.replace(' ','_')}"] = g["X"]
        out_np[f"V_{name.replace(' ','_')}"] = g["V"]
        out_np[f"phi_{name.replace(' ','_')}"] = phi
        out_np[f"edges_{name.replace(' ','_')}"] = np.stack([g["R"], g["C"]])

    # ------- fit on fit region -----------------------------------------
    g = region_graph(fit_name, users, lats, lons, counts, E)
    st = graph_stats(g["n"], g["R"], g["C"])
    phi, bw = spatial_score(g["X"])
    target = np.array([st["ell"], st["C"], st["r"]])
    model = CondModel(g["X"], phi)
    b0 = results[fit_name]["b_direct"]
    x0 = np.array([1.0, 1.0, max(b0, 0.02)])
    log("fitting (lambda, R_km, b) by simulated moments ...")
    th_hat, q_hat = fit_cond(model, target, x0, RNG)
    log(f"fit: lambda={th_hat[0]:.4f} R={th_hat[1]:.4f}km b={th_hat[2]:.4f} "
        f"(psi_graph={th_hat[2]**2:.5f}) q={q_hat:.3e}")
    m_fit = model.sim(*th_hat, RNG, n_rep=60).mean(0)
    log(f"model at fit: ell={m_fit[0]:.3f} C={m_fit[1]:.4f} r={m_fit[2]:.4f} "
        f"vs obs {np.round(target,4)}")
    results["fit"] = dict(theta=[float(x) for x in th_hat],
                          psi_graph=float(th_hat[2] ** 2),
                          q=float(q_hat), model_moments=[float(x) for x in m_fit])

    # ------- ablation: b=0, recalibrate (lambda, R) to (ell, C) --------
    log("ablation: b=0 with (lambda,R) recalibrated to (ell,C) ...")
    th_geo, q_geo = fit_cond(model, target, th_hat, RNG, fix_b=0.0)
    m_geo = model.sim(*th_geo, RNG, n_rep=60).mean(0)
    share = 100.0 * m_geo[2] / target[2]
    log(f"geometry-only: r_geo={m_geo[2]:.4f} vs r_obs={target[2]:.4f} "
        f"-> {share:.1f}% of observed assortativity is geometric baseline")
    results["ablation"] = dict(theta=[float(x) for x in th_geo],
                               model_moments=[float(x) for x in m_geo],
                               geometric_share_pct=float(share))

    # ------- parametric bootstrap CI for psi_graph ----------------------
    log("parametric bootstrap for psi_graph (B=40) ...")
    B = 40
    psis = []
    for bb in range(B):
        mb = model.sim(*th_hat, RNG, n_rep=1)[0]
        thb, _ = fit_cond(model, mb, th_hat, RNG, n_rep=12)
        psis.append(thb[2] ** 2)
        if (bb + 1) % 10 == 0:
            log(f"  bootstrap {bb+1}/{B}")
    psis = np.array(psis)
    results["fit"]["psi_boot_se"] = float(psis.std(ddof=1))
    results["fit"]["psi_boot_q"] = [float(np.quantile(psis, q)) for q in (0.025, 0.975)]
    log(f"psi_graph = {th_hat[2]**2:.5f} +/- {psis.std(ddof=1):.5f} "
        f"(boot 95% {results['fit']['psi_boot_q']})")

    # ------- validation on non-fitted statistics ------------------------
    log("validation: non-fitted statistics on fit region ...")
    def full_sim_stats(theta, n_rep=40):
        lam, Rkm, b = theta
        within = model.pd <= Rkm
        P0, P1 = model.P[:, 0], model.P[:, 1]
        acc = []
        for _ in range(n_rep):
            V = sample_marks_cond(model.phi, b, RNG)
            p = -np.expm1(-lam * V[P0] * V[P1])
            keep = within & (RNG.random(p.size) < p)
            stt = graph_stats(model.n, P0[keep], P1[keep])
            acc.append([stt["deg_var"], stt["H3"] / model.n, stt["H11"] / model.n])
        return np.array(acc)
    obs_extra = np.array([st["deg_var"], st["H3"] / g["n"], st["H11"] / g["n"]])
    sim_fit = full_sim_stats(th_hat)
    sim_geo = full_sim_stats(th_geo)
    results["validation"] = dict(
        observed=dict(deg_var=float(obs_extra[0]), H3n=float(obs_extra[1]), H11n=float(obs_extra[2])),
        model=dict(mean=[float(x) for x in sim_fit.mean(0)], sd=[float(x) for x in sim_fit.std(0)]),
        ablation=dict(mean=[float(x) for x in sim_geo.mean(0)], sd=[float(x) for x in sim_geo.std(0)]))
    log(f"observed extra stats: {np.round(obs_extra,3)}")
    log(f"model:    {np.round(sim_fit.mean(0),3)} +/- {np.round(sim_fit.std(0),3)}")
    log(f"ablation: {np.round(sim_geo.mean(0),3)} +/- {np.round(sim_geo.std(0),3)}")

    # ------- transfer to held-out region --------------------------------
    log("held-out transfer: %s -> %s" % (fit_name, held_name))
    gh = region_graph(held_name, users, lats, lons, counts, E)
    sth = graph_stats(gh["n"], gh["R"], gh["C"])
    phih, _ = spatial_score(gh["X"])
    modelh = CondModel(gh["X"], phih)
    mh = modelh.sim(th_hat[0], th_hat[1], th_hat[2], RNG, n_rep=60).mean(0)
    mh0 = modelh.sim(th_hat[0], th_hat[1], 0.0, RNG, n_rep=60).mean(0)
    results["transfer"] = dict(obs=[sth["ell"], sth["C"], sth["r"]],
                               model=[float(x) for x in mh],
                               ablation=[float(x) for x in mh0])
    log(f"held-out obs: ell={sth['ell']:.3f} C={sth['C']:.4f} r={sth['r']:.4f}")
    log(f"transferred model: {np.round(mh,4)}; ablation: {np.round(mh0,4)}")

    with open(OUT + out_prefix + "_results.json", "w") as fh:
        json.dump(results, fh, indent=1)
    np.savez(OUT + out_prefix + "_arrays.npz", **out_np)
    log("DONE -> " + out_prefix + "_results.json")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2],
         sys.argv[3] if len(sys.argv) > 3 else "brightkite")

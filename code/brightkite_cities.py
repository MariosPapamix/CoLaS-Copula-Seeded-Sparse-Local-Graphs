"""Multi-city Brightkite analysis + baselines + diagnostics.
For each of the 5 prespecified metros: observed moments, Spearman mixing,
degree-preserving configuration null, direct psi, graph-only fit (with
bootstrap), compensated ablation, model-implied channel decomposition
(Cov Lambda + E Gamma over D), and a conditional-design Jacobian
diagnostic at theta-hat.  Plus benchmark information scaling at n=12k.
Writes runs/brightkite_cities.json
"""
import sys, json, time
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
from scipy.stats import rankdata
from brightkite_pipeline import (load_homes_csv, region_graph, graph_stats,
                                 spatial_score, CondModel, fit_cond,
                                 direct_psi, sample_marks_cond, REGIONS)

OUT = _p("runs") + _os.sep
rng = np.random.default_rng(20260813)
t0 = time.time()
def log(m): print("[%6.0fs] %s" % (time.time() - t0, m), flush=True)

users, lats, lons, counts, E = load_homes_csv(
    _p("data", "loc-brightkite_edges.txt.gz"),
    _p("data", "brightkite_homes.csv.gz"))

def spearman_r(deg, R, C):
    rk = rankdata(deg, method="average")
    d1 = np.concatenate([rk[R], rk[C]]); d2 = np.concatenate([rk[C], rk[R]])
    return float(np.corrcoef(d1, d2)[0, 1])

def config_null(n, R, C, n_rep=20, swap_mult=15, rng=None):
    """degree-preserving double-edge swaps; returns mean/sd of (C_hat, r)."""
    base = np.stack([R, C], 1)
    outs = []
    for _ in range(n_rep):
        e = base.copy()
        eset = set((int(a) << 20) | int(b) for a, b in e)
        m = e.shape[0]
        n_sw = swap_mult * m
        idx = rng.integers(0, m, 2 * n_sw)
        done = 0; k = 0
        while done < n_sw and k + 1 < idx.size:
            i, j = idx[k], idx[k + 1]; k += 2
            a, b = e[i]; c, d = e[j]
            if len({a, b, c, d}) < 4:
                continue
            na, nb = (min(a, d), max(a, d)), (min(c, b), max(c, b))
            k1 = (na[0] << 20) | na[1]; k2 = (nb[0] << 20) | nb[1]
            if k1 in eset or k2 in eset:
                continue
            eset.discard((min(a, b) << 20) | max(a, b))
            eset.discard((min(c, d) << 20) | max(c, d))
            eset.add(k1); eset.add(k2)
            e[i] = na; e[j] = nb
            done += 1
        st = graph_stats(n, e[:, 0], e[:, 1])
        outs.append([st["C"], st["r"]])
    outs = np.array(outs)
    return outs.mean(0), outs.std(0)

def model_channels(model, theta, n_rep=30, rng=None):
    """Simulate at theta; compute r decomposition via conditional
    intensities: Lambda_i = sum_j p_ij, Gamma_ij = sum_k p_ik p_jk."""
    lam, Rkm, b = theta
    within = model.pd <= Rkm
    P0, P1 = model.P[within, 0], model.P[within, 1]
    n = model.n
    acc = []
    for _ in range(n_rep):
        V = sample_marks_cond(model.phi, b, rng)
        p = -np.expm1(-lam * V[P0] * V[P1])
        Pm = np.zeros((n, n))
        Pm[P0, P1] = p; Pm[P1, P0] = p
        Lam = Pm.sum(1)
        keep = rng.random(p.size) < p
        R, C = P0[keep], P1[keep]
        if R.size < 10:
            continue
        G2 = Pm @ Pm
        Gam = G2[R, C]
        L1 = np.concatenate([Lam[R], Lam[C]]); L2 = np.concatenate([Lam[C], Lam[R]])
        covL = np.cov(L1, L2)[0, 1]
        D = L1.var() + L1.mean()
        acc.append([covL / D, np.concatenate([Gam, Gam]).mean() / D])
    acc = np.array(acc)
    return acc.mean(0), acc.std(0)

def jac_diag(model, theta, rng, n_rep=40):
    """finite-difference Jacobian of the simulated 3-moment map at theta
    (CRN); returns sigma_min and condition number."""
    def m(th):
        return model.sim(th[0], th[1], max(th[2], 0.0),
                         np.random.default_rng(4242), n_rep).mean(0)
    th = np.array(theta, float)
    J = np.zeros((3, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = np.array([0.05, 0.05, 0.02])[j]
        J[:, j] = (m(th + e) - m(th - e)) / (2 * e[j])
    sv = np.linalg.svd(J, compute_uv=False)
    return dict(sv=[float(x) for x in sv], smin=float(sv[-1]),
                cond=float(sv[0] / max(sv[-1], 1e-12)))

results = {}
for name in REGIONS:
    g = region_graph(name, users, lats, lons, counts, E)
    if g["n"] < 300 or g["R"].size < 200:
        continue
    st = graph_stats(g["n"], g["R"], g["C"])
    deg = np.bincount(g["R"], minlength=g["n"]) + np.bincount(g["C"], minlength=g["n"])
    r_s = spearman_r(deg, g["R"], g["C"])
    (Cn, rn), (Cns, rns) = config_null(g["n"], g["R"], g["C"], rng=rng)
    phi, bw = spatial_score(g["X"])
    b_h, se_b, psi_d, se_pd = direct_psi(g["V"], phi)
    model = CondModel(g["X"], phi)
    target = np.array([st["ell"], st["C"], st["r"]])
    th, q = fit_cond(model, target, np.array([1.0, 1.0, max(b_h, 0.02)]),
                     np.random.default_rng(1))
    th0, q0 = fit_cond(model, target, th, np.random.default_rng(2), fix_b=0.0)
    m_fit = model.sim(*th, np.random.default_rng(3), 60).mean(0)
    m_geo = model.sim(*th0, np.random.default_rng(3), 60).mean(0)
    B = 30
    psis = []
    for bb in range(B):
        mb = model.sim(*th, rng, 1)[0]
        thb, _ = fit_cond(model, mb, th, rng, n_rep=10)
        psis.append(thb[2] ** 2)
    psis = np.array(psis)
    ch_mean, ch_sd = model_channels(model, th, rng=rng)
    jd = jac_diag(model, th, rng)
    results[name] = dict(
        n=st["n"], edges=st["edges"], ell=st["ell"], C=st["C"], r=st["r"],
        r_spearman=r_s, config_null=dict(C=float(Cn), C_sd=float(Cns),
                                         r=float(rn), r_sd=float(rns)),
        b_direct=float(b_h), se_b=float(se_b),
        psi_direct=float(psi_d), se_psi=float(se_pd),
        theta=[float(x) for x in th], psi_graph=float(th[2] ** 2),
        psi_boot_q=[float(np.quantile(psis, x)) for x in (0.025, 0.975)],
        model_moments=[float(x) for x in m_fit],
        r_geo=float(m_geo[2]),
        channels=dict(sort=float(ch_mean[0]), overlap=float(ch_mean[1]),
                      sort_sd=float(ch_sd[0]), overlap_sd=float(ch_sd[1])),
        jac=jd, bw_km=float(bw))
    log(f"{name}: n={st['n']} ell={st['ell']:.2f} C={st['C']:.3f} "
        f"r={st['r']:.3f} r_spear={r_s:.3f} | null C={Cn:.3f} r={rn:.3f} | "
        f"psi_d={psi_d:.4f}({se_pd:.4f}) psi_g={th[2]**2:.4f} "
        f"[{np.quantile(psis,0.025):.3f},{np.quantile(psis,0.975):.3f}] | "
        f"r_geo={m_geo[2]:.3f} ch=({ch_mean[0]:.3f},{ch_mean[1]:.3f}) "
        f"smin={jd['smin']:.4f} cond={jd['cond']:.0f}")

# ---- benchmark information scaling at n=12k ------------------------------
log("benchmark scaling n=12000 ...")
rngd = np.random.default_rng(99)
def metro(lat_c, lon_c, n, cores=3, spread=0.19):
    centers = np.stack([lat_c + rngd.normal(0, 0.08, cores),
                        lon_c + rngd.normal(0, 0.10, cores)], 1)
    w = rngd.dirichlet(np.ones(cores))
    which = rngd.choice(cores, n, p=w)
    return centers[which] + rngd.normal(0, spread * 0.55, (n, 2))
homes12 = metro(37.62, -122.17, 12000)
lat_c = np.radians((37.15 + 38.10) / 2)
X12 = np.stack([(homes12[:, 1] + 122.75) * 111.320 * np.cos(lat_c),
                (homes12[:, 0] - 37.15) * 110.574], 1)
phi12, _ = spatial_score(X12)
V12 = sample_marks_cond(phi12, 0.25, rngd)
from scipy.spatial import cKDTree
tree = cKDTree(X12); Pp = np.array(sorted(tree.query_pairs(2.0)))
pp = -np.expm1(-2.0 * V12[Pp[:, 0]] * V12[Pp[:, 1]])
keep = rngd.random(pp.size) < pp
st12 = graph_stats(12000, Pp[keep, 0], Pp[keep, 1])
model12 = CondModel(X12, phi12, R_max=5.0, max_cand_deg=45.0)
tg12 = np.array([st12["ell"], st12["C"], st12["r"]])
th12, _ = fit_cond(model12, tg12, np.array([1.5, 1.5, 0.2]), np.random.default_rng(7))
ps12 = []
for bb in range(30):
    mb = model12.sim(*th12, rngd, 1)[0]
    thb, _ = fit_cond(model12, mb, th12, rngd, n_rep=8)
    ps12.append(thb[2] ** 2)
ps12 = np.array(ps12)
results["benchmark_n12000"] = dict(
    truth=[2.0, 2.0, 0.0625], theta=[float(x) for x in th12],
    psi_graph=float(th12[2] ** 2),
    psi_boot_q=[float(np.quantile(ps12, x)) for x in (0.025, 0.975)])
log(f"benchmark n=12k: theta={np.round(th12,3)} psi={th12[2]**2:.4f} "
    f"CI=[{np.quantile(ps12,0.025):.4f},{np.quantile(ps12,0.975):.4f}] "
    f"(n=3k CI was [0.001,0.150])")

with open(OUT + "brightkite_cities.json", "w") as fh:
    json.dump(results, fh, indent=1)
log("DONE -> brightkite_cities.json")

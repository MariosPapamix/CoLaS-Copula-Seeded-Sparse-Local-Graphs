"""SF Bay (Brightkite) additional diagnostics:
  A. conditional-degree curves E[D2 | D1=k] (ANND): data vs fitted model
     vs geometry ablation vs configuration null.
  B. fitted latent-distance logistic competitor (estimated likelihood):
     p_ij = sigmoid(alpha - beta * d_ij) over all dyads, MLE, then
     simulated (ell, C, r).
  C. boundary-valid bootstrap specification test: T = n (psi_g - psi_d)^2,
     null distribution by parametric bootstrap at the b=0 conditional fit;
     applied to SF; size check in the certified submodel at psi=0.
Writes runs/sf_extras.json
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
from brightkite_pipeline import (load_homes_csv, region_graph, graph_stats,
                                 spatial_score, CondModel, fit_cond,
                                 direct_psi, sample_marks_cond)
from brightkite_cities_lib import config_null

OUT = _p("runs") + _os.sep
rng = np.random.default_rng(55)
t0 = time.time()
def log(m): print("[%6.0fs] %s" % (time.time() - t0, m), flush=True)

users, lats, lons, counts, E = load_homes_csv(
    _p("data", "loc-brightkite_edges.txt.gz"),
    _p("data", "brightkite_homes.csv.gz"))
g = region_graph("SF Bay", users, lats, lons, counts, E)
st = graph_stats(g["n"], g["R"], g["C"])
phi, _ = spatial_score(g["X"])
n = g["n"]
CT = json.load(open(OUT + "brightkite_cities.json"))
th_fit = CT["SF Bay"]["theta"]; th_geo = [th_fit[0], th_fit[1], 0.0]
res = {}

# ---------- A. ANND curves --------------------------------------------
KMAX = 15
def annd(R, C, n):
    deg = np.bincount(R, minlength=n) + np.bincount(C, minlength=n)
    s = np.zeros(KMAX + 2); c = np.zeros(KMAX + 2)
    d1 = np.concatenate([deg[R], deg[C]]); d2 = np.concatenate([deg[C], deg[R]])
    kk = np.minimum(d1, KMAX + 1)
    np.add.at(s, kk, d2); np.add.at(c, kk, 1)
    with np.errstate(invalid="ignore"):
        return s / np.maximum(c, 1), c

a_obs, c_obs = annd(g["R"], g["C"], n)
model = CondModel(g["X"], phi)
def sim_annd(theta, reps=40):
    lam, Rkm, b = theta
    within = model.pd <= Rkm
    P0, P1 = model.P[within, 0], model.P[within, 1]
    acc = np.zeros((reps, KMAX + 2))
    for rrep in range(reps):
        V = sample_marks_cond(model.phi, b, rng)
        p = -np.expm1(-lam * V[P0] * V[P1])
        keep = rng.random(p.size) < p
        acc[rrep], _ = annd(P0[keep], P1[keep], n)
    return np.nanmean(np.where(acc > 0, acc, np.nan), 0)

a_fit = sim_annd(th_fit)
a_geo = sim_annd(th_geo)
# config-null ANND (single long rewire, 12 reps)
accs = []
base = np.stack([g["R"], g["C"]], 1)
for _ in range(12):
    e = base.copy()
    eset = set((int(a) << 20) | int(b) for a, b in e)
    m = e.shape[0]; idx = rng.integers(0, m, 30 * m)
    done = 0; k = 0
    while done < 15 * m and k + 1 < idx.size:
        i, j = idx[k], idx[k + 1]; k += 2
        a, b = e[i]; c, d = e[j]
        if len({a, b, c, d}) < 4: continue
        na, nb = (min(a, d), max(a, d)), (min(c, b), max(c, b))
        k1 = (na[0] << 20) | na[1]; k2 = (nb[0] << 20) | nb[1]
        if k1 in eset or k2 in eset: continue
        eset.discard((min(a, b) << 20) | max(a, b))
        eset.discard((min(c, d) << 20) | max(c, d))
        eset.add(k1); eset.add(k2)
        e[i] = na; e[j] = nb; done += 1
    aa, _ = annd(e[:, 0], e[:, 1], n)
    accs.append(aa)
a_cfg = np.nanmean(np.where(np.array(accs) > 0, np.array(accs), np.nan), 0)
res["annd"] = dict(k=list(range(KMAX + 2)),
                   obs=[float(x) for x in a_obs],
                   count=[float(x) for x in c_obs],
                   fit=[float(x) for x in a_fit],
                   geo=[float(x) for x in a_geo],
                   cfg=[float(x) for x in a_cfg])
log("A: ANND curves done (k=1..%d+)" % KMAX)

# ---------- B. latent-distance logistic competitor ---------------------
X = g["X"]
iu, ju = np.triu_indices(n, 1)
d = np.linalg.norm(X[iu] - X[ju], axis=1)
y = np.zeros(iu.size, dtype=bool)
key = g["R"].astype(np.int64) * n + g["C"]
kk = iu.astype(np.int64) * n + ju
y[np.isin(kk, key)] = True
ab = np.array([-4.0, 0.05])
for it in range(60):
    eta = ab[0] - ab[1] * d
    p = 1.0 / (1.0 + np.exp(-eta))
    W = p * (1 - p)
    g1 = np.array([(y - p).sum(), -((y - p) * d).sum()])
    H = np.array([[W.sum(), -(W * d).sum()],
                  [-(W * d).sum(), (W * d * d).sum()]])
    step = np.linalg.solve(H, g1)
    ab = ab + step
    if np.max(np.abs(step)) < 1e-10: break
alpha, beta = ab
p_hat = 1.0 / (1.0 + np.exp(-(alpha - beta * d)))
sims = []
for _ in range(40):
    keep = rng.random(p_hat.size) < p_hat
    s2 = graph_stats(n, iu[keep], ju[keep])
    sims.append([s2["ell"], s2["C"], s2["r"]])
sims = np.array(sims)
res["latent_distance"] = dict(alpha=float(alpha), beta=float(beta),
                              moments_mean=[float(x) for x in sims.mean(0)],
                              moments_sd=[float(x) for x in sims.std(0)],
                              observed=[st["ell"], st["C"], st["r"]])
log(f"B: latent-distance fit alpha={alpha:.3f} beta={beta:.4f}/km -> "
    f"model (ell,C,r)={np.round(sims.mean(0),4)} vs obs "
    f"({st['ell']:.2f},{st['C']:.3f},{st['r']:.3f})")

# ---------- C. boundary bootstrap specification test -------------------
# C1: SF application
def sf_pair(theta, rng):
    lam, Rkm, b = theta
    within = model.pd <= Rkm
    P0, P1 = model.P[within, 0], model.P[within, 1]
    V = sample_marks_cond(model.phi, b, rng)
    p = -np.expm1(-lam * V[P0] * V[P1])
    keep = rng.random(p.size) < p
    s2 = graph_stats(n, P0[keep], P1[keep])
    m3 = np.array([s2["ell"], s2["C"], s2["r"]])
    b_d = float(np.mean(np.sqrt(3.0) * (2 * V - 1) * model.phi))
    return m3, b_d

psi_d_obs = CT["SF Bay"]["psi_direct"]
psi_g_obs = CT["SF Bay"]["psi_graph"]
T_obs = n * (psi_g_obs - psi_d_obs) ** 2
B = 200
Tstar = np.empty(B)
for bb in range(B):
    m3s, bds = sf_pair(th_geo, rng)      # null: b = 0 at the conditional fit
    ths, _ = fit_cond(model, m3s, th_geo, rng, n_rep=8)
    Tstar[bb] = n * (max(ths[2], 0.0) ** 2 - bds ** 2) ** 2
    if (bb + 1) % 50 == 0: log(f"C1: bootstrap {bb+1}/{B}")
pval = float(np.mean(Tstar >= T_obs))
res["boundary_test_SF"] = dict(T_obs=float(T_obs), p=pval,
                               crit95=float(np.quantile(Tstar, 0.95)))
log(f"C1: SF boundary bootstrap test: T_obs={T_obs:.5f} p={pval:.3f}")

# C2: size in the certified submodel at psi=0
from colas_sim import Estimator
from colas_moments import MomentMap
from sim_spec_test import gen_full
mm = MomentMap(N=60); est = Estimator(mm)
n_sub, OUTER, INNER = 10000, 200, 99
rej = 0
rng2 = np.random.default_rng(909)
for rep in range(OUTER):
    m3, bd = gen_full(n_sub, 3.0, 9.0, 0.0, rng2)
    thh, _ = est.solve_mom(m3, x0=(2.8, 8.6, 0.02))
    T = n_sub * (max(thh[2], 0.0) - bd ** 2) ** 2
    lam0, eta0 = float(np.clip(thh[0], 1, 5)), float(np.clip(thh[1], 5, 15))
    Ts = np.empty(INNER)
    for ib in range(INNER):
        m3s, bds = gen_full(n_sub, lam0, eta0, 0.0, rng2)
        ths, _ = est.solve_mom(m3s, x0=(lam0, eta0, 0.02))
        Ts[ib] = n_sub * (max(ths[2], 0.0) - bds ** 2) ** 2
    if T > np.quantile(Ts, 0.95): rej += 1
    if (rep + 1) % 40 == 0:
        log(f"C2: size rep {rep+1}/{OUTER} running rate {rej/(rep+1):.3f}")
res["boundary_test_size"] = dict(size=rej / OUTER, outer=OUTER, inner=INNER)
log(f"C2: boundary bootstrap size at psi=0: {rej/OUTER:.3f} (nominal 0.05)")

with open(OUT + "sf_extras.json", "w") as fh:
    json.dump(res, fh, indent=1)
log("DONE -> sf_extras.json")

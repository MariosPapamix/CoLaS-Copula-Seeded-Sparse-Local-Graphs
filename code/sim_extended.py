"""Extended simulation study:
  A. design GRID over the certified box (8 designs incl. weak corners,
     near-boundary, boundary): bias, RMSE, marginal + simultaneous (Wald
     ellipse) coverage, CI length, solver success, channel-share coverage.
  B. five-motif overidentified GMM vs three-moment exactly identified,
     same replicate graphs; J-test (chi^2_2) size.
  C. boundary test size at THREE nuisance configurations.
  D. same-(ell,r), different-mechanism pair: r cannot separate, C does,
     the method recovers the sorting shares.
  E. simultaneous (ellipse) coverage recomputed for the core runs.
Writes runs/sim_r3.npz
"""
import sys, time
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
from scipy import stats
from colas_sim import simulate, sim_targets, Estimator
from colas_moments import MomentMap, channels, channel_ci, five_moment_map

OUT = _p("runs") + _os.sep
mm = MomentMap(N=60)
est = Estimator(mm)
res = {}
t0 = time.time()
def log(m): print("[%6.0fs] %s" % (time.time() - t0, m), flush=True)

CHI3_95 = stats.chi2.ppf(0.95, 3)
CRIT5 = 2.705543

def sim5(n, th, rng):
    s = simulate(n, *th, rng)
    y5 = np.array([s["ell"], s["T"]/n, s["S2"]/n, s["H3"]/n, s["H11"]/n])
    return np.array([s["ell"], s["C"], s["r"]]), y5

# ---------------- A. grid --------------------------------------------------
DESIGNS = [(1.5, 6.0, 0.02), (1.5, 14.0, 0.12), (3.0, 9.0, 0.0835),
           (4.5, 6.0, 0.05), (4.5, 14.0, 0.14), (2.0, 12.0, 0.005),
           (3.0, 9.0, 0.15), (2.5, 7.0, 0.04)]
n = 10000
grid_rows = []
for d_i, th0 in enumerate(DESIGNS):
    rng = np.random.default_rng(5000 + d_i)
    B = 500
    Mb = np.array([sim_targets(n, *th0, rng) for _ in range(B)])
    Sigma = np.cov(Mb.T) * n
    REP = 300
    TH = np.empty((REP, 3)); OK = np.zeros(REP, bool)
    COV = np.zeros((REP, 3), bool); ELL = np.zeros(REP, bool)
    LEN = np.empty(REP); SHARE_COV = np.zeros(REP, bool)
    rs0, ro0 = channels(mm, *th0)
    share0 = rs0 / (rs0 + ro0)
    for rep in range(REP):
        m = sim_targets(n, *th0, rng)
        th, ok = est.solve_mom(m, x0=(max(th0[0]*0.8,0.5), th0[1]*0.9, max(th0[2],0.01)))
        TH[rep] = th; OK[rep] = ok
        ci, V = est.ci_delta(th, Sigma, n)
        COV[rep] = np.abs(th - np.array(th0)) <= ci
        LEN[rep] = 2 * ci[2]
        dv = th - np.array(th0)
        try:
            W = n * dv @ np.linalg.solve(V * n, dv)
            ELL[rep] = W <= CHI3_95
        except np.linalg.LinAlgError:
            ELL[rep] = False
        (rs, ro, sh), se = channel_ci(mm, th, V)
        SHARE_COV[rep] = abs(sh - share0) <= 1.96 * se[2]
    bias = TH.mean(0) - np.array(th0)
    rmse = np.sqrt(((TH - np.array(th0))**2).mean(0))
    grid_rows.append(dict(theta=th0, bias=bias, rmse=rmse,
                          cov_marg=COV.mean(0), cov_ell=ELL.mean(),
                          len_psi=np.median(LEN), ok=OK.mean(),
                          share0=share0, share_cov=SHARE_COV.mean()))
    log(f"A: design {th0}: cov={COV.mean(0).round(3)} ell={ELL.mean():.3f} "
        f"rmse_psi={rmse[2]:.4f} share0={share0:.3f} share_cov={SHARE_COV.mean():.3f} ok={OK.mean():.3f}")
res["grid"] = grid_rows

# ---------------- B. five-moment GMM vs three-moment ----------------------
th0 = (3.0, 9.0, 0.0835)
rng = np.random.default_rng(777001)
B = 1200
Y5b = np.empty((B, 5))
for b in range(B):
    _, Y5b[b] = sim5(n, th0, rng)
Om5 = np.cov(Y5b.T) * n
W5 = np.linalg.inv(Om5)
M3b = np.array([sim_targets(n, *th0, np.random.default_rng(888000 + b)) for b in range(500)])
Sig3 = np.cov(M3b.T) * n
log("B: banks done")

def gmm5(y5, x0):
    x = np.array(x0, float)
    lo = np.array([0.3, 2.0, -0.10]); hi = np.array([6.5, 25.0, 0.163])
    def mu(th): return five_moment_map(mm, *th)
    for _ in range(60):
        m0 = mu(x); d = y5 - m0
        G = np.zeros((5, 3))
        for j in range(3):
            e = np.zeros(3); e[j] = 1e-5 * max(1, abs(x[j]))
            G[:, j] = (mu(x + e) - mu(x - e)) / (2 * e[j])
        g = G.T @ W5 @ d
        H = G.T @ W5 @ G
        try:
            dx = np.linalg.solve(H, g)
        except np.linalg.LinAlgError:
            return x, np.inf, None
        xn = np.clip(x + dx, lo, hi)
        if np.linalg.norm(xn - x) < 1e-10:
            x = xn; break
        x = xn
    d = y5 - mu(x)
    Q = float(d @ W5 @ d)
    G = np.zeros((5, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = 1e-5 * max(1, abs(x[j]))
        G[:, j] = (mu(x + e) - mu(x - e)) / (2 * e[j])
    V = np.linalg.inv(G.T @ W5 @ G) / n
    return x, Q, V

REP = 400
TH3 = np.empty((REP, 3)); TH5 = np.empty((REP, 3))
COV5 = np.zeros((REP, 3), bool); JSTAT = np.empty(REP)
rng = np.random.default_rng(31415)
for rep in range(REP):
    m3, y5 = sim5(n, th0, rng)
    th3, _ = est.solve_mom(m3)
    th5, Q, V5 = gmm5(y5, th3)
    TH3[rep] = th3; TH5[rep] = th5
    JSTAT[rep] = n * Q
    if V5 is not None:
        se = np.sqrt(np.maximum(np.diag(V5), 0))
        COV5[rep] = np.abs(th5 - np.array(th0)) <= 1.96 * se
res["fm_th3"] = TH3; res["fm_th5"] = TH5; res["fm_J"] = JSTAT
res["fm_cov5"] = COV5
rm3 = np.sqrt(((TH3 - np.array(th0))**2).mean(0))
rm5 = np.sqrt(((TH5 - np.array(th0))**2).mean(0))
log(f"B: RMSE 3-moment {rm3.round(4)} vs 5-moment {rm5.round(4)}; "
    f"5m coverage {COV5.mean(0).round(3)}; J-test size at 5%: "
    f"{np.mean(JSTAT > stats.chi2.ppf(0.95, 2)):.3f}")

# ---------------- C. boundary size at three nuisances ---------------------
bnd = {}
for (l0, e0) in ((2.0, 7.0), (3.0, 9.0), (4.0, 12.0)):
    rng = np.random.default_rng(int(l0 * 1000 + e0))
    Mb = np.array([sim_targets(n, l0, e0, 0.0, rng) for _ in range(400)])
    Sg = np.cov(Mb.T) * n; Wn = np.linalg.inv(Sg)
    REP = 600; D = np.empty(REP)
    for rep in range(REP):
        m = sim_targets(n, l0, e0, 0.0, rng)
        th, ok = est.solve_mom(m, x0=(l0 * 0.95, e0 * 0.95, 0.02))
        if th[2] <= 0:
            D[rep] = 0.0
        else:
            _, q0 = est.solve_null(m, Wn, x0=(th[0], th[1]))
            D[rep] = n * q0
    bnd[(l0, e0)] = dict(size=float(np.mean(D > CRIT5)), mass0=float(np.mean(D <= 1e-12)))
    log(f"C: nuisance ({l0},{e0}): size={bnd[(l0,e0)]['size']:.3f} mass0={bnd[(l0,e0)]['mass0']:.3f}")
res["boundary_nuisance"] = {str(k): v for k, v in bnd.items()}

# ---------------- D. same-(ell,r) different-mechanism pair -----------------
thA = np.array([2.2, 11.0, 0.13])
tA = mm.targets(*thA)
target_lr = np.array([float(tA["ell"]), float(tA["r"])])
x = np.array([3.0, 9.0]); psiB = 0.005
for _ in range(60):
    t = mm.targets(x[0], x[1], psiB)
    F = np.array([float(t["ell"]), float(t["r"])]) - target_lr
    if np.max(np.abs(F)) < 1e-11: break
    J = t["J"].reshape(3, 3)[[0, 2]][:, :2]
    x = x - np.linalg.solve(J, F)
thB = np.array([x[0], x[1], psiB])
tB = mm.targets(*thB)
rsA, roA = channels(mm, *thA); rsB, roB = channels(mm, *thB)
log(f"D: pair matched on (ell,r)=({target_lr[0]:.3f},{target_lr[1]:.3f}): "
    f"A={thA.round(3)} C={float(tA['C']):.4f} share={rsA/(rsA+roA):.3f} | "
    f"B={thB.round(3)} C={float(tB['C']):.4f} share={rsB/(rsB+roB):.3f}")
REP = 400
outAB = {}
for tag, th in (("A", thA), ("B", thB)):
    rng = np.random.default_rng(hash(tag) % 2**31)
    M = np.array([sim_targets(n, *th, rng) for _ in range(REP)])
    TH = np.empty((REP, 3)); SH = np.empty(REP)
    for rep in range(REP):
        thh, ok = est.solve_mom(M[rep], x0=(2.5, 9.0, 0.05))
        TH[rep] = thh
        rs, ro = channels(mm, *np.clip(thh, [0.5, 3, 0.0], [6, 20, 0.16]))
        SH[rep] = rs / (rs + ro)
    outAB[tag] = dict(m=M, th=TH, share=SH, theta=th,
                      C_true=float(mm.targets(*th)["C"]),
                      share_true=(rsA/(rsA+roA) if tag == "A" else rsB/(rsB+roB)))
    log(f"D: {tag}: rhat={M[:,2].mean():.4f}(sd {M[:,2].std():.4f}) "
        f"Chat={M[:,1].mean():.4f}(sd {M[:,1].std():.4f}) share_hat={SH.mean():.3f}(sd {SH.std():.3f})")
res["pairA_m"] = outAB["A"]["m"]; res["pairB_m"] = outAB["B"]["m"]
res["pairA_share"] = outAB["A"]["share"]; res["pairB_share"] = outAB["B"]["share"]
res["pairA_theta"] = thA; res["pairB_theta"] = thB
res["pair_shares_true"] = np.array([outAB["A"]["share_true"], outAB["B"]["share_true"]])
res["pair_C_true"] = np.array([outAB["A"]["C_true"], outAB["B"]["C_true"]])

# ---------------- E. simultaneous coverage for core runs -------------------
R2 = np.load(OUT + "sim_results.npz", allow_pickle=True)
sim_cov = {}
for nn in (1000, 10000, 100000):
    TH = R2[f"cov_theta_{nn}"]; Sg = R2[f"bank_Sigma_{nn}"]
    ELL = []
    for th in TH:
        _, J = est.targets(th)
        Ji = np.linalg.inv(J)
        V = Ji @ Sg @ Ji.T / nn
        dv = th - np.array([3.0, 9.0, 0.0835])
        try:
            ELL.append(dv @ np.linalg.solve(V, dv) <= CHI3_95)
        except np.linalg.LinAlgError:
            ELL.append(False)
    sim_cov[nn] = float(np.mean(ELL))
    log(f"E: n={nn}: Wald-ellipse simultaneous coverage {sim_cov[nn]:.3f} "
        f"(naive all-3-marginal joint {R2[f'cov_hit_{nn}'].all(1).mean():.3f})")
res["ellipse_cov"] = np.array([sim_cov[1000], sim_cov[10000], sim_cov[100000]])
res["joint_naive"] = np.array([R2[f"cov_hit_{nn}"].all(1).mean() for nn in (1000, 10000, 100000)])

np.savez(OUT + "sim_r3.npz", **{k: v for k, v in res.items() if not isinstance(v, (dict, list))})
import json
with open(OUT + "sim_r3_extra.json", "w") as fh:
    json.dump({"grid": [{k: (v.tolist() if isinstance(v, np.ndarray) else v)
                         for k, v in row.items()} for row in grid_rows],
               "boundary_nuisance": res["boundary_nuisance"]}, fh, indent=1)
log("DONE -> sim_r3.npz / sim_r3_extra.json")

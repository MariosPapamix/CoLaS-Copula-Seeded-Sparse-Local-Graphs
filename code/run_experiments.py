"""Full simulation study for the paper.
Stages:
  A. oracle covariance banks at theta0 for n in {1e3,1e4,1e5}
  B. coverage/bias of MoM/GMM with oracle-Sigma CIs (+ empirical bias of moments)
  C. feasible arm at n=1e4 (own pilot + R=150 simulated covariance per rep)
  D. boundary psi=0: null distribution of D_n, n^{1/4} amplitude rate
  E. power curve at n=1e4 over psi grid
  F. misspecification: triangular/bump kernels, two-harmonic marks
Writes runs/sim_results.npz  (+ progress log lines)
"""
import sys, time, json
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
from colas_sim import simulate, sim_targets, Estimator
from colas_moments import MomentMap

OUT = _p("runs") + _os.sep
mm = MomentMap(N=60, cache=OUT + "coeffs_N60.pkl")
est = Estimator(mm)

TH0 = (3.0, 9.0, 0.0835)
NULL = (3.0, 9.0, 0.0)
CRIT5 = 2.705543           # P(0.5 chi2_0 + 0.5 chi2_1 > c) = 0.05
results = {}
t_start = time.time()
def log(msg):
    print("[%7.0fs] %s" % (time.time() - t_start, msg), flush=True)

# ---------------- A: oracle banks --------------------------------------
banks = {}
for n, B in ((1000, 2000), (10000, 2000), (100000, 400)):
    rng = np.random.default_rng(1234 + n)
    M = np.empty((B, 3))
    for b in range(B):
        M[b] = sim_targets(n, *TH0, rng)
    banks[n] = M
    Sigma = np.cov(M.T) * n
    results[f"bank_mean_{n}"] = M.mean(0)
    results[f"bank_Sigma_{n}"] = Sigma
    log(f"A: bank n={n} B={B} done; mean={np.round(M.mean(0),5)}")
np.savez(OUT + "sim_banks.npz", **{f"M_{n}": banks[n] for n in banks})

# finite-n bias of moments (bank mean - theory), scaled by n
th_theory, _ = est.targets(TH0)
for n in banks:
    bias = banks[n].mean(0) - th_theory
    results[f"moment_bias_{n}"] = bias
    log(f"A: n={n} moment bias = {np.round(bias,6)} ; n*bias = {np.round(n*bias,3)}")

# ---------------- B: coverage / bias -----------------------------------
cover = {}
for n, REP in ((1000, 400), (10000, 400), (100000, 100)):
    rng = np.random.default_rng(777 + n)
    Sigma = results[f"bank_Sigma_{n}"]
    TH = np.empty((REP, 3)); OK = np.zeros(REP, bool); COV = np.zeros((REP, 3), bool)
    SES = np.empty((REP, 3))
    for rep in range(REP):
        m = sim_targets(n, *TH0, rng)
        th, ok = est.solve_mom(m)
        TH[rep] = th; OK[rep] = ok
        ci, _ = est.ci_delta(th, Sigma, n)
        SES[rep] = ci / 1.96
        COV[rep] = np.abs(th - np.array(TH0)) <= ci
    results[f"cov_theta_{n}"] = TH
    results[f"cov_hit_{n}"] = COV
    results[f"cov_ok_{n}"] = OK
    results[f"cov_se_{n}"] = SES
    log(f"B: n={n} coverage lam/eta/psi = {COV.mean(0).round(3)}, "
        f"bias = {np.round(TH.mean(0)-np.array(TH0),5)}, ok={OK.mean():.3f}")

# ---------------- C: feasible arm at n=1e4 ------------------------------
n, REP, RPI = 10000, 100, 150
rng = np.random.default_rng(31337)
TH = np.empty((REP, 3)); COV = np.zeros((REP, 3), bool)
for rep in range(REP):
    m = sim_targets(n, *TH0, rng)
    th, ok = est.solve_mom(m)
    Mp = np.empty((RPI, 3))
    for b in range(RPI):
        Mp[b] = sim_targets(n, th[0], th[1], max(th[2], 0.0), rng)
    Sig_p = np.cov(Mp.T) * n
    ci, _ = est.ci_delta(th, Sig_p, n)
    TH[rep] = th
    COV[rep] = np.abs(th - np.array(TH0)) <= ci
results["feas_theta"] = TH; results["feas_hit"] = COV
log(f"C: feasible n=1e4 coverage = {COV.mean(0).round(3)}")

# ---------------- D: boundary null -------------------------------------
for n, REP in ((1000, 2000), (10000, 1000), (100000, 200)):
    rng = np.random.default_rng(999 + n)
    Sigma = results[f"bank_Sigma_{n}"] if n in (1000, 10000, 100000) else None
    # null banks: covariance at NULL, not TH0 -> build quick bank
    Bn = 600 if n < 100000 else 200
    Mn = np.empty((Bn, 3))
    for b in range(Bn):
        Mn[b] = sim_targets(n, *NULL, rng)
    Sig0 = np.cov(Mn.T) * n
    W = np.linalg.inv(Sig0)
    results[f"null_Sigma_{n}"] = Sig0
    D = np.empty(REP); PSI = np.empty(REP)
    for rep in range(REP):
        m = sim_targets(n, *NULL, rng)
        th, ok = est.solve_mom(m, x0=(2.8, 8.5, 0.02))
        PSI[rep] = th[2]
        if th[2] <= 0.0:
            D[rep] = 0.0
        else:
            _, q0 = est.solve_null(m, W, x0=(th[0], th[1]))
            D[rep] = n * q0
    results[f"null_D_{n}"] = D
    results[f"null_psi_{n}"] = PSI
    log(f"D: null n={n}: P(D=0)={np.mean(D<=1e-12):.3f} "
        f"(theory 0.5), P(D>crit)={np.mean(D>CRIT5):.4f} (nominal 0.05)")

# I_eff at the null (for the rate law constant)
_, J0 = est.targets(NULL)
Sig0 = results["null_Sigma_10000"]
Om_inv = np.linalg.inv(Sig0)
Ga = J0[:, :2]; gp = J0[:, 2]
P = Ga @ np.linalg.inv(Ga.T @ Om_inv @ Ga) @ Ga.T @ Om_inv
rpsi = gp - P @ gp
results["I_eff_null"] = float(rpsi @ Om_inv @ rpsi)
log(f"D: I_eff at null = {results['I_eff_null']:.4f}")

# ---------------- E: power at n=1e4 ------------------------------------
n, REP = 10000, 500
psis = [0.0, 0.005, 0.01, 0.02, 0.04, 0.0835]
Sig0 = results["null_Sigma_10000"]; W = np.linalg.inv(Sig0)
POW = []
rng = np.random.default_rng(2468)
for ps in psis:
    rej = 0
    for rep in range(REP):
        m = sim_targets(n, 3.0, 9.0, ps, rng)
        th, ok = est.solve_mom(m, x0=(2.8, 8.5, max(ps, 0.02)))
        if th[2] <= 0:
            continue
        _, q0 = est.solve_null(m, W, x0=(th[0], th[1]))
        if n * q0 > CRIT5:
            rej += 1
    POW.append(rej / REP)
    log(f"E: power psi={ps}: {rej/REP:.3f}")
results["power_psis"] = np.array(psis)
results["power"] = np.array(POW)

# ---------------- F: misspecification ----------------------------------
n, REP = 10000, 300
for tag, kw in (("tri", dict(kernel="tri")),
                ("bump", dict(kernel="bump")),
                ("twoharm", dict(harmonics=[(np.sqrt(0.0835/2), 1, 0.0),
                                            (np.sqrt(0.0835/2), 2, 1.0)]))):
    rng = np.random.default_rng(hash(tag) % 2**31)
    TH = np.empty((REP, 3))
    for rep in range(REP):
        m = sim_targets(n, 3.0, 9.0, 0.0835 if tag != "twoharm" else 0.0, rng, **kw)
        th, ok = est.solve_mom(m)
        TH[rep] = th
    results[f"miss_{tag}"] = TH
    log(f"F: {tag}: mean theta = {TH.mean(0).round(4)}  (true psi_tot=0.0835)")

np.savez(OUT + "sim_results.npz", **results)
log("ALL DONE -> sim_results.npz")

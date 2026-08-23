"""Hausman graph-vs-direct specification test (size + power),
heavier robustness (long-range ties, heavy-tailed mark transform),
collapsed-displacement baseline constant.  Writes runs/sim_r4.npz"""
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
from colas_sim import sample_marks, Estimator, motif_stats
from colas_moments import MomentMap

OUT = _p("runs") + _os.sep
mm = MomentMap(N=60)
est = Estimator(mm)
TH0 = (3.0, 9.0, 0.0835)
n = 10000
t0 = time.time()
def log(m): print("[%6.0fs] %s" % (time.time() - t0, m), flush=True)

def gen_full(n, lam, eta, psi, rng, mark_transform=None, frac_random=0.0):
    """graph stats + direct estimate; optional misspecifications."""
    X = rng.random(n)
    order = np.argsort(X); xs = X[order]
    V = sample_marks(xs, psi, rng)
    W = mark_transform(V) if mark_transform else V
    w = eta / n
    idx = np.arange(n)
    hi = np.searchsorted(xs, xs + w, side="right")
    cnt = hi - idx - 1
    tot = int(cnt.sum())
    I = np.repeat(idx, cnt)
    start = np.concatenate(([0], np.cumsum(cnt)[:-1]))
    J = I + 1 + (np.arange(tot) - np.repeat(start, cnt))
    over = xs + w - 1.0
    m2 = over > 0
    hi2 = np.where(m2, np.searchsorted(xs, np.maximum(over, 0.0), side="right"), 0)
    hi2 = np.minimum(hi2, idx)
    tot2 = int(hi2.sum())
    if tot2:
        I2 = np.repeat(idx, hi2)
        st2 = np.concatenate(([0], np.cumsum(hi2)[:-1]))
        J2 = np.arange(tot2) - np.repeat(st2, hi2)
        I = np.concatenate([I, I2]); J = np.concatenate([J, J2])
    p = -np.expm1(-lam * W[I] * W[J])
    keep = rng.random(p.size) < p
    R, C = I[keep], J[keep]
    if frac_random > 0:
        m_extra = int(frac_random * R.size)
        A = rng.integers(0, n, m_extra); B = rng.integers(0, n, m_extra)
        ok = A != B
        R = np.concatenate([R, np.minimum(A[ok], B[ok])])
        C = np.concatenate([C, np.maximum(A[ok], B[ok])])
    deg = np.bincount(R, minlength=n) + np.bincount(C, minlength=n)
    s = motif_stats(R, C, deg, n)
    b_dir = float(np.mean(np.sqrt(3.0) * (2 * V - 1) * np.sqrt(2.0) * np.cos(2 * np.pi * xs)))
    return np.array([s["ell"], s["C"], s["r"]]), b_dir

def psi_pair(m3, b_dir):
    th, ok = est.solve_mom(m3)
    return max(th[2], 0.0), b_dir ** 2

# oracle bank for sigma_Delta and the individual variances
B = 800
rng = np.random.default_rng(424242)
diffs = np.empty(B); pg = np.empty(B); pd_ = np.empty(B)
for b in range(B):
    m3, bd = gen_full(n, *TH0, rng)
    g, d = psi_pair(m3, bd)
    pg[b] = g; pd_[b] = d; diffs[b] = g - d
var_delta = float(np.var(diffs))
log(f"bank: sd(psi_g)={pg.std():.4f} sd(psi_d)={pd_.std():.4f} "
    f"sd(diff)={diffs.std():.4f} corr={np.corrcoef(pg,pd_)[0,1]:.3f}")

CHI1_95 = stats.chi2.ppf(0.95, 1)
def run_arm(tag, REP=500, **kw):
    rng = np.random.default_rng(hash(tag) % 2**31)
    T = np.empty(REP); PG = np.empty(REP); PD = np.empty(REP)
    for rep in range(REP):
        m3, bd = gen_full(n, *TH0, rng, **kw)
        g, d = psi_pair(m3, bd)
        PG[rep] = g; PD[rep] = d
        T[rep] = (g - d) ** 2 / var_delta
    rej = float(np.mean(T > CHI1_95))
    log(f"{tag}: reject={rej:.3f}  psi_g={PG.mean():.4f}({PG.std():.4f}) "
        f"psi_d={PD.mean():.4f}")
    return dict(rej=rej, psi_g=PG, psi_d=PD, T=T)

res = {}
res["size"] = run_arm("H0 correct spec")
res["lr2"] = run_arm("2% long-range ties", frac_random=0.02)
res["lr5"] = run_arm("5% long-range ties", frac_random=0.05)
res["heavy"] = run_arm("heavy-tail mark transform",
                       mark_transform=lambda v: 0.1 + 1.4 * v * v)

np.savez(OUT + "sim_r4.npz",
         var_delta=var_delta,
         size_rej=res["size"]["rej"], lr2_rej=res["lr2"]["rej"],
         lr5_rej=res["lr5"]["rej"], heavy_rej=res["heavy"]["rej"],
         size_T=res["size"]["T"],
         lr2_psig=res["lr2"]["psi_g"], lr5_psig=res["lr5"]["psi_g"],
         heavy_psig=res["heavy"]["psi_g"],
         size_psig=res["size"]["psi_g"], size_psid=res["size"]["psi_d"])
log("DONE -> sim_r4.npz")

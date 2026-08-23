"""
colas_sim.py -- fast simulator + estimator for the CoLaS d=1 hard-range
submodel (rho=1, eps_n = 1/n, W=V, one harmonic m=1, phi(x)=sqrt2 cos 2pi x).

Graph law: X_i ~ U[0,1); V_i | X_i ~ density 1 + sqrt(psi) a(v) phi(x),
a(v)=sqrt3(2v-1); edge for i<j independent with prob
   1 - exp(-lambda V_i V_j k(n*d_T(X_i,X_j)))          [rho_n = 1 exactly]
with k = 1{|z|<=eta} by default (misspecification kernels optional).

Statistics: E_n, T_n, S2_n, H3, H11, mean degree ell=2E/n, transitivity
C=3T/S2, Newman endpoint-degree correlation r over directed edges.
"""
from __future__ import annotations
import numpy as np
import scipy.sparse as sp

# ------------------------- mark sampling ------------------------------
SQRT3 = np.sqrt(3.0)
SQRT2 = np.sqrt(2.0)

def sample_marks(X, psi, rng, harmonics=None):
    """V | X for one-harmonic (default) or a list [(b, m, phase)]."""
    if harmonics is None:
        s = SQRT3 * np.sqrt(psi) * SQRT2 * np.cos(2 * np.pi * X)
    else:
        f = np.zeros_like(X)
        for (b, m, ph) in harmonics:
            f += b * SQRT2 * np.cos(2 * np.pi * m * X - ph)
        s = SQRT3 * f
    u = rng.random(X.shape)
    small = np.abs(s) < 1e-9
    s_safe = np.where(small, 1.0, s)
    disc = (1.0 - s_safe) ** 2 + 4.0 * s_safe * u
    v = (-(1.0 - s_safe) + np.sqrt(np.maximum(disc, 0.0))) / (2.0 * s_safe)
    return np.where(small, u, np.clip(v, 0.0, 1.0))

# ------------------------- graph generation ---------------------------
def gen_graph(n, lam, eta, psi, rng, kernel="hard", harmonics=None):
    """Returns (rows, cols, deg) with rows<cols realized edges."""
    X = rng.random(n)
    order = np.argsort(X)
    xs = X[order]
    V = sample_marks(xs, psi, rng, harmonics)
    w = eta / n
    idx = np.arange(n)
    # right neighbors within w (linear part)
    hi = np.searchsorted(xs, xs + w, side="right")
    cnt = hi - idx - 1
    tot = int(cnt.sum())
    I = np.repeat(idx, cnt)
    start = np.concatenate(([0], np.cumsum(cnt)[:-1]))
    J = I + 1 + (np.arange(tot) - np.repeat(start, cnt))
    # wraparound pairs
    over = xs + w - 1.0
    m2 = over > 0
    hi2 = np.where(m2, np.searchsorted(xs, np.maximum(over, 0.0), side="right"), 0)
    hi2 = np.minimum(hi2, idx)          # j < i always (top block vs bottom block)
    cnt2 = hi2
    tot2 = int(cnt2.sum())
    if tot2:
        I2 = np.repeat(idx, cnt2)
        start2 = np.concatenate(([0], np.cumsum(cnt2)[:-1]))
        J2 = np.arange(tot2) - np.repeat(start2, cnt2)
        I = np.concatenate([I, I2]); J = np.concatenate([J, J2])
    # distances and kernel
    d = np.abs(xs[J] - xs[I])
    d = np.minimum(d, 1.0 - d) * n          # tangent-scale distance
    if kernel == "hard":
        kv = 1.0
        keep0 = d <= eta
    elif kernel == "tri":
        kv = 2.0 * np.maximum(1.0 - d / eta, 0.0)
        keep0 = d <= eta
    elif kernel == "bump":
        kv = np.exp(-2.0 * d * d / (eta * eta)) / 0.5981620  # normalized: int = 2*eta
        keep0 = d <= eta
    else:
        raise ValueError(kernel)
    p = -np.expm1(-lam * V[I] * V[J] * kv)
    keep = keep0 & (rng.random(p.shape) < p)
    R, C = I[keep], J[keep]
    deg = np.bincount(R, minlength=n) + np.bincount(C, minlength=n)
    return R, C, deg, n

def motif_stats(R, C, deg, n):
    E_n = R.size
    S2 = 0.5 * float((deg.astype(np.float64) * (deg - 1)).sum())
    H3 = float((deg.astype(np.float64) ** 3).sum())
    A = sp.csr_matrix((np.ones(2 * R.size, dtype=np.float32),
                       (np.concatenate([R, C]), np.concatenate([C, R]))), shape=(n, n))
    # triangles (blocked to bound memory)
    T = 0.0
    step = max(1, n // 8) if n > 30000 else n
    for s in range(0, n, step):
        blk = A[s:s + step]
        T += float((blk @ A).multiply(blk).sum())
    T /= 6.0
    Ad = A @ deg.astype(np.float64)
    H11 = float(deg.astype(np.float64) @ Ad)
    # Newman r over directed edges
    d1 = np.concatenate([deg[R], deg[C]]).astype(np.float64)
    d2 = np.concatenate([deg[C], deg[R]]).astype(np.float64)
    if d1.size >= 2 and d1.std() > 0:
        r = float(np.corrcoef(d1, d2)[0, 1])
    else:
        r = 0.0
    ell = 2.0 * E_n / n
    Ch = 3.0 * T / S2 if S2 > 0 else 0.0
    return dict(E=E_n, T=T, S2=S2, H3=H3, H11=H11, ell=ell, C=Ch, r=r)

def simulate(n, lam, eta, psi, rng, kernel="hard", harmonics=None):
    R, C, deg, n = gen_graph(n, lam, eta, psi, rng, kernel, harmonics)
    return motif_stats(R, C, deg, n)

def sim_targets(n, lam, eta, psi, rng, **kw):
    s = simulate(n, lam, eta, psi, rng, **kw)
    return np.array([s["ell"], s["C"], s["r"]])

# ------------------------- estimation ---------------------------------
class Estimator:
    def __init__(self, mm):
        self.mm = mm    # MomentMap

    def targets(self, th):
        t = self.mm.targets(th[0], th[1], th[2])
        return np.array([float(t["ell"]), float(t["C"]), float(t["r"])]), t["J"].reshape(3, 3)

    def solve_mom(self, m_obs, x0=(2.5, 8.0, 0.05), itmax=80):
        x = np.array(x0, float)
        lo = np.array([0.3, 2.0, -0.10]); hi = np.array([6.5, 25.0, 0.163])
        F = None
        for _ in range(itmax):
            f, J = self.targets(x)
            F = f - m_obs
            if np.max(np.abs(F)) < 1e-11:
                return x, True
            try:
                dx = np.linalg.solve(J, F)
            except np.linalg.LinAlgError:
                return x, False
            step = 1.0
            for _ in range(30):
                xn = np.clip(x - step * dx, lo, hi)
                fn, _ = self.targets(xn)
                if np.linalg.norm(fn - m_obs) < np.linalg.norm(F) or step < 1e-6:
                    break
                step *= 0.5
            x = xn
        return x, bool(np.max(np.abs(F)) < 1e-7)

    def solve_null(self, m_obs, W, x0=(2.5, 8.0), itmax=100):
        """min over (lam, eta) of (m - M3(lam,eta,0))' W (m - M3(lam,eta,0))."""
        x = np.array(x0, float)
        lo = np.array([0.3, 2.0]); hi = np.array([6.5, 25.0])
        def qval(x2):
            f, _ = self.targets((x2[0], x2[1], 0.0))
            d = f - m_obs
            return float(d @ W @ d), d, f
        q, d, f = qval(x)
        for _ in range(itmax):
            _, J = self.targets((x[0], x[1], 0.0))
            G = J[:, :2]
            g = 2.0 * G.T @ W @ d
            H = 2.0 * G.T @ W @ G
            try:
                dx = np.linalg.solve(H, g)
            except np.linalg.LinAlgError:
                break
            step = 1.0
            improved = False
            for _ in range(30):
                xn = np.clip(x - step * dx, lo, hi)
                qn, dn, fn = qval(xn)
                if qn < q - 1e-16:
                    x, q, d, f = xn, qn, dn, fn
                    improved = True
                    break
                step *= 0.5
            if not improved or np.linalg.norm(g) < 1e-13:
                break
        return x, q

    def ci_delta(self, th_hat, Sigma_m, n, level=1.96):
        _, J = self.targets(th_hat)
        Ji = np.linalg.inv(J)
        V = Ji @ Sigma_m @ Ji.T / n
        se = np.sqrt(np.maximum(np.diag(V), 0.0))
        return se * level, V

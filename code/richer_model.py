"""Richer conditional class: two-scale connection kernel +
node-level degree propensities + one-pass triadic closure, fitted by
simulated moment matching to eight targets
    (ell, C, r, q25, q50, q75 of edge length, degree CV, p99 degree)
with the GOF-before-attribution discipline of the main pipeline.

Model, conditional on observed home coordinates X and density score phi:
  marks    V_i ~ h_b(. | phi_i)                    (same sorting family)
  propensities s_i = exp(sigma Z_i - sigma^2/2),   Z_i iid N(0,1)
  kernel   kappa(d) = 1{d <= R} + w exp(-d / L)
  base     p_ij = 1 - exp(-lam s_i s_j V_i V_j kappa(d_ij))
  closure  one pass: non-adjacent pair with m common neighbors is added
           with probability 1 - (1 - tau)^m
Parameters theta = (lam, R, w, L, sigma, tau, b); b is SIGNED (the
conditional design can carry linear-in-b information).
"""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy.sparse as sp
from scipy.spatial.distance import pdist
from brightkite_pipeline import sample_marks_cond

TARGET_KEYS = ["ell", "C", "r", "q25", "q50", "q75", "deg_cv", "p99"]

class RichModel:
    def __init__(self, X, phi):
        self.X = np.asarray(X, float)
        self.phi = np.asarray(phi, float)
        self.n = self.X.shape[0]
        n = self.n
        self.D = pdist(self.X).astype(np.float32)      # condensed pair dists
        iu = np.triu_indices(n, 1)
        self.I = iu[0].astype(np.int32)
        self.J = iu[1].astype(np.int32)

    def _stats(self, R_idx, C_idx, d_edge):
        n = self.n
        deg = (np.bincount(R_idx, minlength=n)
               + np.bincount(C_idx, minlength=n)).astype(float)
        E_n = R_idx.size
        if E_n < 10:
            return None
        S2 = 0.5 * float((deg * (deg - 1)).sum())
        A = sp.csr_matrix((np.ones(2 * E_n, dtype=np.float32),
                           (np.concatenate([R_idx, C_idx]),
                            np.concatenate([C_idx, R_idx]))), shape=(n, n))
        T = float((A @ A).multiply(A).sum()) / 6.0
        d1 = np.concatenate([deg[R_idx], deg[C_idx]])
        d2 = np.concatenate([deg[C_idx], deg[R_idx]])
        r = float(np.corrcoef(d1, d2)[0, 1]) if d1.std() > 0 else 0.0
        q = np.percentile(d_edge, [25, 50, 75]) if d_edge.size else [0, 0, 0]
        return dict(ell=2.0 * E_n / n, C=3.0 * T / S2 if S2 > 0 else 0.0,
                    r=r, q25=float(q[0]), q50=float(q[1]), q75=float(q[2]),
                    deg_cv=float(deg.std() / max(deg.mean(), 1e-9)),
                    p99=float(np.percentile(deg, 99)),
                    deg=deg, A=A, edges=E_n)

    def sim(self, theta, rng, n_rep=1, full=False):
        lam, R, w, L, sigma, tau, b = theta
        out = []
        for _ in range(n_rep):
            V = sample_marks_cond(self.phi, b, rng)
            s = np.exp(sigma * rng.standard_normal(self.n) - 0.5 * sigma ** 2)
            wV = (s * V).astype(np.float32)
            kap = (self.D <= R).astype(np.float32)
            kap += np.float32(w) * np.exp(-self.D / np.float32(max(L, 1e-3)))
            rate = np.float32(lam) * wV[self.I] * wV[self.J] * kap
            p = -np.expm1(-rate)
            keep = rng.random(p.size) < p
            R_idx = self.I[keep].astype(np.int64)
            C_idx = self.J[keep].astype(np.int64)
            d_edge = self.D[keep].astype(float)
            # one-pass closure
            if tau > 1e-9 and R_idx.size:
                n = self.n
                A = sp.csr_matrix((np.ones(2 * R_idx.size, dtype=np.float32),
                                   (np.concatenate([R_idx, C_idx]),
                                    np.concatenate([C_idx, R_idx]))),
                                  shape=(n, n))
                M = sp.triu(A @ A, k=1).tocoo()
                if M.nnz:
                    have = set((int(a) << 21) | int(b)
                               for a, b in zip(R_idx, C_idx))
                    cand = [(a, b, m) for a, b, m in
                            zip(M.row, M.col, M.data)
                            if ((int(a) << 21) | int(b)) not in have]
                    if cand:
                        ca = np.array([c[0] for c in cand], dtype=np.int64)
                        cb = np.array([c[1] for c in cand], dtype=np.int64)
                        cm = np.array([c[2] for c in cand], dtype=float)
                        pc = 1.0 - (1.0 - tau) ** cm
                        add = rng.random(pc.size) < pc
                        if add.any():
                            na, nb = ca[add], cb[add]
                            R_idx = np.concatenate([R_idx, na])
                            C_idx = np.concatenate([C_idx, nb])
                            dd = np.linalg.norm(self.X[na] - self.X[nb],
                                                axis=1)
                            d_edge = np.concatenate([d_edge, dd])
            st = self._stats(R_idx, C_idx, d_edge)
            if st is None:
                continue
            if full:
                out.append(st)
            else:
                out.append([st[k] for k in TARGET_KEYS])
        return out if full else np.array(out)

def rel_resid(m, target):
    """Standardized residual vector: relative for scale quantities,
    absolute/0.05 for r (its natural scale)."""
    v = []
    for k, key in enumerate(TARGET_KEYS):
        t, x = target[k], m[k]
        if key == "r":
            v.append((x - t) / 0.03)
        elif key == "C":
            v.append((x - t) / max(abs(t), 1e-6) / 0.05)
        else:
            v.append((x - t) / max(abs(t), 1e-6) / 0.10)   # 10% units
    return np.array(v)

def fit_rich(model, target, rng_seed=999, sweeps=4, n_rep=6, x0=None,
             fix_b=None, verbose=True):
    """Hierarchical calibration + Nelder-Mead polish, CRN throughout."""
    lam, R, w, L, sigma, tau, b = x0 if x0 is not None else \
        (2.0, 1.0, 0.03, 18.0, 1.0, 0.04, 0.0)
    if fix_b is not None:
        b = fix_b

    def mom(theta, reps=n_rep, seed=rng_seed):
        m = model.sim(theta, np.random.default_rng(seed), reps)
        return m.mean(0) if len(m) else np.full(len(TARGET_KEYS), np.nan)

    idx = {k: i for i, k in enumerate(TARGET_KEYS)}
    for sw in range(sweeps):
        # 1. lam -> ell
        for _ in range(8):
            m = mom((lam, R, w, L, sigma, tau, b))
            f = target[idx["ell"]] / max(m[idx["ell"]], 1e-9)
            if abs(f - 1) < 0.01:
                break
            lam = float(np.clip(lam * f, 1e-3, 500.0))
        # 2. (R, w, L) -> quantiles
        for _ in range(6):
            m = mom((lam, R, w, L, sigma, tau, b))
            fR = target[idx["q25"]] / max(m[idx["q25"]], 1e-9)
            fL = target[idx["q75"]] / max(m[idx["q75"]], 1e-9)
            fw = target[idx["q50"]] / max(m[idx["q50"]], 1e-9)
            R = float(np.clip(R * fR ** 0.6, 0.05, 45.0))
            L = float(np.clip(L * fL ** 0.6, 0.5, 200.0))
            # more mid-length mass -> raise w (q50 low means too local)
            w = float(np.clip(w * fw ** 0.8, 1e-4, 5.0))
        # 3. sigma -> deg_cv
        for _ in range(6):
            m = mom((lam, R, w, L, sigma, tau, b))
            f = target[idx["deg_cv"]] / max(m[idx["deg_cv"]], 1e-9)
            if abs(f - 1) < 0.02:
                break
            sigma = float(np.clip(sigma * f ** 0.8, 0.05, 3.0))
        # 4. tau -> C
        for _ in range(6):
            m = mom((lam, R, w, L, sigma, tau, b))
            f = target[idx["C"]] / max(m[idx["C"]], 1e-9)
            if abs(f - 1) < 0.02:
                break
            tau = float(np.clip(tau * f ** 0.9, 1e-4, 0.9))
        # 5. b -> r (signed, only if not fixed)
        if fix_b is None:
            for _ in range(5):
                m = mom((lam, R, w, L, sigma, tau, b))
                gap = target[idx["r"]] - m[idx["r"]]
                if abs(gap) < 0.004:
                    break
                b = float(np.clip(b + 1.5 * gap, -0.55, 0.55))
        if verbose:
            m = mom((lam, R, w, L, sigma, tau, b))
            rr = rel_resid(m, target)
            print(f"[rich sweep {sw}] theta=({lam:.3f},{R:.3f},{w:.4f},"
                  f"{L:.1f},{sigma:.3f},{tau:.4f},{b:+.3f}) "
                  f"maxres={np.abs(rr).max():.2f}", flush=True)

    # Nelder-Mead polish on all free parameters
    from scipy.optimize import minimize
    free = [0, 1, 2, 3, 4, 5] + ([] if fix_b is not None else [6])
    th0 = np.array([lam, R, w, L, sigma, tau, b])

    def obj(z):
        th = th0.copy()
        th[free] = z
        th[0] = abs(th[0]); th[1] = np.clip(th[1], 0.05, 45)
        th[2] = np.clip(th[2], 1e-4, 5); th[3] = np.clip(th[3], 0.5, 200)
        th[4] = np.clip(th[4], 0.05, 3); th[5] = np.clip(th[5], 1e-4, 0.9)
        th[6] = np.clip(th[6], -0.55, 0.55)
        m = mom(tuple(th), reps=n_rep)
        if np.isnan(m).any():
            return 1e6
        return float((rel_resid(m, target) ** 2).sum())

    res = minimize(obj, th0[free], method="Nelder-Mead",
                   options=dict(maxiter=int(os.environ.get("RICH_NM", "160")), xatol=1e-3, fatol=1e-3))
    th = th0.copy(); th[free] = res.x
    th[0] = abs(th[0]); th[1] = np.clip(th[1], 0.05, 45)
    th[2] = np.clip(th[2], 1e-4, 5); th[3] = np.clip(th[3], 0.5, 200)
    th[4] = np.clip(th[4], 0.05, 3); th[5] = np.clip(th[5], 1e-4, 0.9)
    th[6] = np.clip(th[6], -0.55, 0.55)
    if verbose:
        m = mom(tuple(th), reps=max(n_rep, 12))
        print(f"[rich polish] theta={np.round(th,4).tolist()} "
              f"obj={res.fun:.3f} maxres={np.abs(rel_resid(m,target)).max():.2f}",
              flush=True)
    return th

def gof(model, theta, target, B=100, seed=20260816):
    """Bootstrap GOF at theta: mean, sd, z per target."""
    rng = np.random.default_rng(seed)
    ms = model.sim(tuple(theta), rng, B)
    mean, sd = ms.mean(0), ms.std(0)
    z = (np.array(target) - mean) / np.maximum(sd, 1e-9)
    return dict(mean=mean.tolist(), sd=sd.tolist(), z=z.tolist(),
                max_abs_z=float(np.abs(z).max()), B=B)

"""
colas_moments.py
================
Exact polynomial reduction of the CoLaS hard-range submodel moment map
(d=1, rho=1, W=V~Unif(0,1), one harmonic m=1, hard window kernel).

The five mark-motif integrals over the unit cube,

    E = int H12 w2,           (edge)
    S = int H12 H13 w3,       (wedge, center 1)
    T = int H12 H13 H23 w3,   (triangle)
    P = int H12 H13 H24 w4,   (3-path 3-1-2-4)
    Q = int H12 H13 H14 w4,   (3-star, center 1)

with H_ij = 1 - exp(-lambda u_i u_j) and w_q the location-averaged
one-harmonic mark density, reduce exactly (after truncating
1-e^{-x} = P_N(x) + R_N(x)) to polynomials

    M(lambda, psi) = sum_s lambda^s (a0[s] + a1[s] psi + a2[s] psi^2),

using the monomial moments p_d = 1/(d+1),
qb_d = d/((d+1)(d+2))  ( = (1/sqrt3) int u^d a(u) du ),
and

    Omega_{v,d}(psi) = prod_i p_{d_i}
                     + 3 psi sum_{i<j} qb_{d_i} qb_{d_j} prod_{k not in {i,j}} p_{d_k}
                     + 1{v=4} (27/2) psi^2 prod_i qb_{d_i},

where d_i = sum of edge exponents at vertex i.  Signs: each term of the
expanded product of P_N's carries (-1)^{sum n_e + m} lambda^{sum n_e} / prod n_e!.

Truncation error (from SI S3, eqs. product-value/lambda/psi-tail):
    |prod H - prod P_N|          <= m eps0 (1+eps0)^{m-1}
    |d/dlam (prod H - prod P_N)| <= m {eps1 (1+eps0)^{m-1} + (m-1) eps0 (1+eps0)^{(m-2)+}}
    integrated psi-derivative err <= B_v m eps0 (1+eps0)^{m-1},  (B2,B3,B4)=(3,9,23)
with eps0 = L^{N+1}/(N+1)!, eps1 = L^N/N!  (L = max lambda), and the value
errors multiplied by 2^v because 0 <= w_v <= 2^v.

Observable map (tangent scale):
    ell = 2 eta E
    C   = (3/4) T / S                       (free of eta)
    r   = num/den,  num = 2 eta (P E - S^2) + (3/4) T E
                    den = 2 eta (Q E - S^2) + S E
Reduced map at fixed ell (2 eta = ell / E):
    rt  = num'/den', num' = ell (P E - S^2) + (3/4) T E^2
                     den' = ell (Q E - S^2) + S E^2

This module provides float (numpy) evaluation of all values and first
partials; exact Fraction coefficients are exported for the certified
verifier (colas_certificate_v2.py).
"""
from __future__ import annotations

import numpy as np
import os as _os
_BASE = _os.environ.get("COLAS_BASE") or _os.path.dirname(
    _os.path.dirname(_os.path.abspath(__file__)))
def _p(*parts):
    q = _os.path.join(_BASE, *parts)
    _os.makedirs(_os.path.dirname(q), exist_ok=True) if parts[:1] in (("runs",), ("figures",), ("revision",)) else None
    return q

from fractions import Fraction
from math import factorial

# ----------------------------------------------------------------------
# Motif definitions: edge lists on vertices 0..v-1
# ----------------------------------------------------------------------
MOTIFS = {
    "E": dict(v=2, edges=[(0, 1)]),
    "S": dict(v=3, edges=[(0, 1), (0, 2)]),
    "T": dict(v=3, edges=[(0, 1), (0, 2), (1, 2)]),
    "P": dict(v=4, edges=[(0, 1), (0, 2), (1, 3)]),
    "Q": dict(v=4, edges=[(0, 1), (0, 2), (0, 3)]),
}

def _p(d):   return Fraction(1, d + 1)
def _qb(d):  return Fraction(d, (d + 1) * (d + 2))


def build_coeffs(N: int, motif: str):
    """Exact Fraction coefficient arrays (a0, a1, a2), index = total lambda degree s,
    for truncation order N per edge.  a2 nonzero only for v=4 motifs."""
    spec = MOTIFS[motif]
    v, edges = spec["v"], spec["edges"]
    m = len(edges)
    smax = N * m
    a0 = [Fraction(0)] * (smax + 1)
    a1 = [Fraction(0)] * (smax + 1)
    a2 = [Fraction(0)] * (smax + 1)

    # iterate over exponent tuples
    import itertools
    for tup in itertools.product(range(1, N + 1), repeat=m):
        s = sum(tup)
        coef = Fraction((-1) ** (s + m))
        for n_e in tup:
            coef /= factorial(n_e)
        d = [0] * v
        for (i, j), n_e in zip(edges, tup):
            d[i] += n_e
            d[j] += n_e
        ps = [_p(di) for di in d]
        qs = [_qb(di) for di in d]
        prod_p = Fraction(1)
        for x in ps:
            prod_p *= x
        # psi^1 term: 3 * sum_{i<j} qb_i qb_j * prod_{k not in {i,j}} p_k
        t1 = Fraction(0)
        for i in range(v):
            for j in range(i + 1, v):
                term = qs[i] * qs[j]
                for k in range(v):
                    if k != i and k != j:
                        term *= ps[k]
                t1 += term
        t1 *= 3
        a0[s] += coef * prod_p
        a1[s] += coef * t1
        if v == 4:
            prod_q = Fraction(1)
            for x in qs:
                prod_q *= x
            a2[s] += coef * Fraction(27, 2) * prod_q
    return a0, a1, a2


class MomentMap:
    """Float evaluator for (E,S,T,P,Q), the observable map (ell, C, r),
    the reduced map rt(lambda,psi,ell), and all first partials."""

    def __init__(self, N: int = 60, cache: str | None = None):
        self.N = N
        if cache is None:
            import os
            cache = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 "coeffs_N%d.pkl" % N)
        if cache is not None:
            import os, pickle
            if os.path.exists(cache):
                with open(cache, "rb") as fh:
                    self.coef = pickle.load(fh)
                if self.coef.get("N") == N:
                    self._to_float()
                    return
        self.coef = {"N": N}
        for k in MOTIFS:
            self.coef[k] = build_coeffs(N, k)
        if cache is not None:
            import pickle
            with open(cache, "wb") as fh:
                pickle.dump(self.coef, fh)
        self._to_float()

    def _to_float(self):
        self.fc = {}
        for k in MOTIFS:
            a0, a1, a2 = self.coef[k]
            self.fc[k] = (np.array([float(x) for x in a0]),
                          np.array([float(x) for x in a1]),
                          np.array([float(x) for x in a2]))

    # -- raw motifs ----------------------------------------------------
    def motif(self, k, lam, psi, deriv=False):
        """Return M (and if deriv: dM/dlam, dM/dpsi); lam, psi broadcastable arrays."""
        a0, a1, a2 = self.fc[k]
        s = np.arange(len(a0))
        lam = np.asarray(lam, dtype=float)
        psi = np.asarray(psi, dtype=float)
        L = lam[..., None] ** s          # (..., smax+1)
        c = a0 + psi[..., None] * a1 + (psi[..., None] ** 2) * a2
        val = (L * c).sum(-1)
        if not deriv:
            return val
        dLam = (np.where(s > 0, s, 0) * lam[..., None] ** np.maximum(s - 1, 0) * c).sum(-1)
        dPsi = (L * (a1 + 2 * psi[..., None] * a2)).sum(-1)
        return val, dLam, dPsi

    def all_motifs(self, lam, psi):
        out = {}
        for k in MOTIFS:
            v, dl, dp = self.motif(k, lam, psi, deriv=True)
            out[k] = (v, dl, dp)
        return out

    # -- observable map and Jacobian ----------------------------------
    def targets(self, lam, eta, psi):
        """Return dict with ell, C, r and full 3x3 Jacobian wrt (lam, eta, psi)."""
        lam = np.asarray(lam, float); eta = np.asarray(eta, float); psi = np.asarray(psi, float)
        M = self.all_motifs(lam, psi)
        E, El, Ep = M["E"]; S, Sl, Sp = M["S"]; T, Tl, Tp = M["T"]
        P, Pl, Pp = M["P"]; Q, Ql, Qp = M["Q"]

        ell = 2 * eta * E
        ell_l = 2 * eta * El; ell_e = 2 * E; ell_p = 2 * eta * Ep

        C = 0.75 * T / S
        C_l = 0.75 * (Tl * S - T * Sl) / S**2
        C_p = 0.75 * (Tp * S - T * Sp) / S**2

        num = 2 * eta * (P * E - S**2) + 0.75 * T * E
        den = 2 * eta * (Q * E - S**2) + S * E
        num_l = 2 * eta * (Pl * E + P * El - 2 * S * Sl) + 0.75 * (Tl * E + T * El)
        num_p = 2 * eta * (Pp * E + P * Ep - 2 * S * Sp) + 0.75 * (Tp * E + T * Ep)
        num_e = 2 * (P * E - S**2)
        den_l = 2 * eta * (Ql * E + Q * El - 2 * S * Sl) + Sl * E + S * El
        den_p = 2 * eta * (Qp * E + Q * Ep - 2 * S * Sp) + Sp * E + S * Ep
        den_e = 2 * (Q * E - S**2)
        r = num / den
        r_l = (num_l * den - num * den_l) / den**2
        r_e = (num_e * den - num * den_e) / den**2
        r_p = (num_p * den - num * den_p) / den**2

        J = np.stack([np.stack([ell_l, ell_e, ell_p], -1),
                      np.stack([C_l, np.zeros_like(C_l), C_p], -1),
                      np.stack([r_l, r_e, r_p], -1)], -2)
        return dict(ell=ell, C=C, r=r, J=J,
                    motifs=dict(E=E, S=S, T=T, P=P, Q=Q),
                    den=den)

    # -- reduced map at fixed ell (global-injectivity object) ---------
    def reduced(self, lam, psi, ell):
        """rt(lam,psi;ell) = num'/den' and the gate quantities
        C_lambda, Delta = rt_psi*C_l - rt_lam*C_p, profile slope Delta/C_l."""
        lam = np.asarray(lam, float); psi = np.asarray(psi, float); ell = np.asarray(ell, float)
        M = self.all_motifs(lam, psi)
        E, El, Ep = M["E"]; S, Sl, Sp = M["S"]; T, Tl, Tp = M["T"]
        P, Pl, Pp = M["P"]; Q, Ql, Qp = M["Q"]

        C_l = 0.75 * (Tl * S - T * Sl) / S**2
        C_p = 0.75 * (Tp * S - T * Sp) / S**2

        A = P * E - S**2; B = Q * E - S**2
        A_l = Pl * E + P * El - 2 * S * Sl; A_p = Pp * E + P * Ep - 2 * S * Sp
        B_l = Ql * E + Q * El - 2 * S * Sl; B_p = Qp * E + Q * Ep - 2 * S * Sp
        num = ell * A + 0.75 * T * E**2
        den = ell * B + S * E**2
        num_l = ell * A_l + 0.75 * (Tl * E**2 + 2 * T * E * El)
        num_p = ell * A_p + 0.75 * (Tp * E**2 + 2 * T * E * Ep)
        den_l = ell * B_l + Sl * E**2 + 2 * S * E * El
        den_p = ell * B_p + Sp * E**2 + 2 * S * E * Ep
        rt = num / den
        rt_l = (num_l * den - num * den_l) / den**2
        rt_p = (num_p * den - num * den_p) / den**2
        Delta = rt_p * C_l - rt_l * C_p
        return dict(rt=rt, C_l=C_l, C_p=C_p, rt_l=rt_l, rt_p=rt_p,
                    Delta=Delta, slope=Delta / C_l, den=den)

    # -- inverse (Newton) ---------------------------------------------
    def solve(self, ell_obs, C_obs, r_obs, x0=(3.0, 9.0, 0.08), tol=1e-12, itmax=60):
        """Solve M3(lam,eta,psi) = (ell,C,r).  Returns (theta, ok)."""
        x = np.array(x0, float)
        target = np.array([ell_obs, C_obs, r_obs], float)
        for _ in range(itmax):
            t = self.targets(x[0], x[1], x[2])
            F = np.array([t["ell"], t["C"], t["r"]], float) - target
            if np.max(np.abs(F)) < tol:
                return x, True
            J = t["J"].reshape(3, 3)
            try:
                dx = np.linalg.solve(J, F)
            except np.linalg.LinAlgError:
                return x, False
            step = 1.0
            # damped Newton with positivity constraints (psi may go <=0 freely: allow, clamp later)
            for _ in range(40):
                xn = x - step * dx
                if xn[0] > 1e-6 and xn[1] > 1e-6:
                    break
                step *= 0.5
            x = x - step * dx
        return x, np.max(np.abs(F)) < 1e-8


if __name__ == "__main__":
    mm = MomentMap(N=60, cache=_p("runs", "coeffs_N60.pkl"))
    # cross-check against the frozen v1 certificate at the K_star box
    t = mm.targets(3.0, 9.0, 0.0835)
    J = t["J"].reshape(3, 3)
    detJ = np.linalg.det(J)
    A2 = np.array([[J[0, 0], J[0, 1]], [J[1, 0], J[1, 1]]])
    detA = np.linalg.det(A2)
    smin = np.linalg.svd(J, compute_uv=False)[-1]
    red = mm.reduced(3.0, 0.0835, t["ell"])
    print("center of K_star (3, 9, 0.0835):")
    print("  ell=%.6f  C=%.6f  r=%.6f" % (t["ell"], t["C"], t["r"]))
    print("  detJ=%.6f (v1 frozen range [-0.078436,-0.055595])" % detJ)
    print("  ||J||_F^2=%.5f (v1 sup <= 6.29724)" % (J * J).sum())
    print("  sigma_min=%.6f (v1 lower bd 0.017657)" % smin)
    print("  2E=ell_eta=%.6f (v1 min >= 0.881935)" % J[0, 1])
    print("  C_lambda=%.6f (v1 min >= 0.061628)" % J[1, 0])
    print("  profile slope detJ/detA=%.6f (v1 min >= 0.999809)" % (detJ / detA))
    print("  reduced-map slope Delta/C_l=%.6f (must match)" % red["slope"])
    print("  E=%.6f S=%.6f T=%.6f P=%.6f Q=%.6f" % tuple(t["motifs"][k] for k in "ESTPQ"))


# ----------------------------------------------------------------------
# Channel decomposition (added in revision round 3):
#   r = r_sort + r_overlap with
#   r_sort    = 2 eta (P E - S^2) / D      (endpoint-intensity covariance)
#   r_overlap = (3/4) T E / D              (shared-neighbor covariance)
#   D         = 2 eta (Q E - S^2) + S E
# ----------------------------------------------------------------------
def channels(mm, lam, eta, psi):
    M = mm.all_motifs(np.asarray(lam, float), np.asarray(psi, float))
    E = M["E"][0]; S = M["S"][0]; T = M["T"][0]; P = M["P"][0]; Q = M["Q"][0]
    D = 2 * eta * (Q * E - S**2) + S * E
    r_sort = 2 * eta * (P * E - S**2) / D
    r_ov = 0.75 * T * E / D
    return float(r_sort), float(r_ov)

def channel_ci(mm, th_hat, V_theta, eps=1e-5):
    """Delta-method SEs for (r_sort, r_overlap, share=r_sort/r) at th_hat,
    given V_theta = Cov(theta_hat)."""
    th_hat = np.asarray(th_hat, float)
    def f(th):
        rs, ro = channels(mm, *th)
        return np.array([rs, ro, rs / (rs + ro)])
    base = f(th_hat)
    G = np.zeros((3, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = eps * max(1.0, abs(th_hat[j]))
        G[:, j] = (f(th_hat + e) - f(th_hat - e)) / (2 * e[j])
    V = G @ V_theta @ G.T
    se = np.sqrt(np.maximum(np.diag(V), 0.0))
    return base, se

def five_moment_map(mm, lam, eta, psi):
    """Closed-form limit means (ell, t, s2, h3, h11)."""
    M = mm.all_motifs(np.asarray(lam, float), np.asarray(psi, float))
    E = M["E"][0]; S = M["S"][0]; T = M["T"][0]; P = M["P"][0]; Q = M["Q"][0]
    ell = 2 * eta * E
    t = 0.5 * eta**2 * T
    s2 = 2 * eta**2 * S
    h3 = 8 * eta**3 * Q + 12 * eta**2 * S + 2 * eta * E
    h11 = 8 * eta**3 * P + 3 * eta**2 * T + 8 * eta**2 * S + 2 * eta * E
    return np.array([float(ell), float(t), float(s2), float(h3), float(h11)])

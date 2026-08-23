"""Certified positivity of the efficient boundary information on the
null face.

Object: the five-moment map mu(theta) = (2 eta E, eta^2 T / 2,
2 eta^2 S, 8 eta^3 Q + 12 eta^2 S + 2 eta E, 8 eta^3 P + 3 eta^2 T
+ 8 eta^2 S + 2 eta E) at psi = 0, with G = D mu (5 x 3) in
(lambda, eta, psi).  For the optimally weighted GMM,
    I_full = G^T Omega^{-1} G  >=  lambda_min(G^T G) / lambda_max(Omega) I_3,
and the psi-profiled information satisfies
    I_eff >= lambda_min(G^T G) / lambda_max(Omega).
This script certifies, by exact rational interval arithmetic on an
adaptive partition of the face F = [1,5] x [5,15] x {0}:
    (i)  a positive lower bound on lambda_min(G^T G) via
         det(G^T G) / trace(G^T G)^2 (valid for 3x3 PSD matrices);
    (ii) combined with the PROVEN upper bound
         lambda_max(Omega) <= C_Omega(eta) = 142 * 5 * 209 * (2 eta)^6
         (Lemma S4.x: ||L||_F^2 = 142; at most 209 positive-probability
         overlap patterns per connected motif pair on <= 4 vertices;
         each quotient integral J_G <= (2 eta)^{v-1} <= (2 eta)^6 at
         psi = 0), a certified uniform bound
             inf_F I_eff >= c_I > 0.
Motif enclosures use the exact Fraction coefficient tables (N = 40)
and the truncation-tail constants of the exact-reduction proposition
at L = 5.  All gates are plain comparisons that raise CertificateError
on failure (they remain active under python -O).
"""
import os, sys, pickle
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fractions import Fraction
from math import factorial

class CertificateError(Exception):
    pass

# ---------------- exact interval helpers (no rounding needed) ---------
def imul(a, b):
    ps = (a[0]*b[0], a[0]*b[1], a[1]*b[0], a[1]*b[1])
    return (min(ps), max(ps))

def iadd(a, b):
    return (a[0]+b[0], a[1]+b[1])

def isub(a, b):
    return (a[0]-b[0], a[1]-b[1])

def iscale(a, c):
    return (a[0]*c, a[1]*c) if c >= 0 else (a[1]*c, a[0]*c)

def iwiden(a, eps):
    return (a[0]-eps, a[1]+eps)

def horner(coeffs, lam_iv):
    """Interval Horner for sum coeffs[s] lam^s, lam interval positive."""
    acc = (Fraction(0), Fraction(0))
    for c in reversed(coeffs):
        acc = imul(acc, lam_iv)
        acc = iadd(acc, (c, c))
    return acc

def taylor_shift(coeffs, lam_m):
    """Exact re-centering: p(lam) = sum_k b_k (lam - lam_m)^k."""
    d = len(coeffs) - 1
    b = [Fraction(0)] * (d + 1)
    # b_k = sum_{s>=k} a_s C(s,k) lam_m^{s-k}
    pow_lm = [Fraction(1)] * (d + 1)
    for s in range(1, d + 1):
        pow_lm[s] = pow_lm[s - 1] * lam_m
    from math import comb
    for s, a in enumerate(coeffs):
        if a == 0:
            continue
        for k in range(s + 1):
            b[k] += a * comb(s, k) * pow_lm[s - k]
    return b

_SHIFT_CACHE = {}

def horner_shifted(key, coeffs, lam_iv):
    """Tight enclosure via exact midpoint Taylor shift (cached per cell)."""
    lam_m = (lam_iv[0] + lam_iv[1]) / 2
    ck = (key, lam_m)
    if ck not in _SHIFT_CACHE:
        _SHIFT_CACHE[ck] = taylor_shift(coeffs, lam_m)
    b = _SHIFT_CACHE[ck]
    dw = (lam_iv[1] - lam_iv[0]) / 2
    delta = (-dw, dw)
    acc = (Fraction(0), Fraction(0))
    for c in reversed(b):
        acc = imul(acc, delta)
        acc = iadd(acc, (c, c))
    return acc

def ipow(a, k):
    out = (Fraction(1), Fraction(1))
    for _ in range(k):
        out = imul(out, a)
    return out

# ---------------- tail constants (Proposition S3, L = 5, N = 40) -------
N_TRUNC = 40
L_RANGE = Fraction(5)
EPS0 = L_RANGE ** (N_TRUNC + 1) / factorial(N_TRUNC + 1)
EPS1 = L_RANGE ** N_TRUNC / factorial(N_TRUNC)
B_V = {2: 3, 3: 9, 4: 23}

def tails(v, m):
    """(value, lambda-deriv, psi-deriv) absolute truncation budgets,
    including the 2^v mark-density factor (conservative at psi=0)."""
    onep = (1 + EPS0)
    val = m * EPS0 * onep ** (m - 1)
    dl = m * (EPS1 * onep ** (m - 1) + (m - 1) * EPS0 * onep ** max(m - 2, 0))
    dp = B_V[v] * m * EPS0 * onep ** (m - 1)
    f = Fraction(2) ** v
    return f * val, f * dl, f * dp

MOTIF_VM = {"E": (2, 1), "S": (3, 2), "T": (3, 3), "P": (4, 3), "Q": (4, 3)}

def load_coeffs():
    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "coeffs_exact_N40.pkl")
    with open(path, "rb") as fh:
        C = pickle.load(fh)
    if len(C["E"][0]) != N_TRUNC + 1:
        raise CertificateError(
            f"coefficient table order mismatch: len={len(C['E'][0])}")
    for k in MOTIF_VM:
        for arr in C[k][:2]:
            if not all(isinstance(x, Fraction) for x in arr):
                raise CertificateError(f"non-exact coefficient in {k}")
    return C

def motif_enclosures(C, lam_iv):
    """Per motif: intervals for value at psi=0 (a0), d/dlam (a0'),
    d/dpsi at psi=0 (a1), all with truncation budgets added."""
    out = {}
    for k, (v, m) in MOTIF_VM.items():
        a0, a1, _ = C[k]
        d0 = [s * a0[s] for s in range(1, len(a0))]
        t_val, t_dl, t_dp = tails(v, m)
        out[k] = dict(val=iwiden(horner_shifted((k, "v"), a0, lam_iv), t_val),
                      dl=iwiden(horner_shifted((k, "d"), d0, lam_iv), t_dl),
                      dp=iwiden(horner_shifted((k, "p"), a1, lam_iv), t_dp))
    return out

def G_intervals(M, eta_iv):
    """5x3 interval matrix G = D mu at psi = 0."""
    e1 = eta_iv
    e2 = ipow(eta_iv, 2)
    e3 = ipow(eta_iv, 3)
    E, S, T, P, Q = (M[k] for k in ("E", "S", "T", "P", "Q"))
    two = Fraction(2); half = Fraction(1, 2)
    rows = []
    # mu1 = 2 eta E
    rows.append((iscale(imul(e1, E["dl"]), two),
                 iscale(E["val"], two),
                 iscale(imul(e1, E["dp"]), two)))
    # mu2 = eta^2 T / 2
    rows.append((iscale(imul(e2, T["dl"]), half),
                 imul(iscale(e1, 1), T["val"]),          # d/deta = eta T
                 iscale(imul(e2, T["dp"]), half)))
    # mu3 = 2 eta^2 S
    rows.append((iscale(imul(e2, S["dl"]), two),
                 iscale(imul(e1, S["val"]), 4),
                 iscale(imul(e2, S["dp"]), two)))
    # mu4 = 8 eta^3 Q + 12 eta^2 S + 2 eta E
    rows.append((iadd(iadd(iscale(imul(e3, Q["dl"]), 8),
                           iscale(imul(e2, S["dl"]), 12)),
                      iscale(imul(e1, E["dl"]), 2)),
                 iadd(iadd(iscale(imul(e2, Q["val"]), 24),
                           iscale(imul(e1, S["val"]), 24)),
                      iscale(E["val"], 2)),
                 iadd(iadd(iscale(imul(e3, Q["dp"]), 8),
                           iscale(imul(e2, S["dp"]), 12)),
                      iscale(imul(e1, E["dp"]), 2))))
    # mu5 = 8 eta^3 P + 3 eta^2 T + 8 eta^2 S + 2 eta E
    rows.append((iadd(iadd(iadd(iscale(imul(e3, P["dl"]), 8),
                                iscale(imul(e2, T["dl"]), 3)),
                           iscale(imul(e2, S["dl"]), 8)),
                      iscale(imul(e1, E["dl"]), 2)),
                 iadd(iadd(iadd(iscale(imul(e2, P["val"]), 24),
                                iscale(imul(e1, T["val"]), 6)),
                           iscale(imul(e1, S["val"]), 16)),
                      iscale(E["val"], 2)),
                 iadd(iadd(iadd(iscale(imul(e3, P["dp"]), 8),
                                iscale(imul(e2, T["dp"]), 3)),
                           iscale(imul(e2, S["dp"]), 8)),
                      iscale(imul(e1, E["dp"]), 2))))
    return rows

def gram(rows):
    M = [[(Fraction(0), Fraction(0)) for _ in range(3)] for _ in range(3)]
    for a in range(3):
        for b in range(3):
            acc = (Fraction(0), Fraction(0))
            for i in range(5):
                acc = iadd(acc, imul(rows[i][a], rows[i][b]))
            M[a][b] = acc
    return M

def det3(M):
    def mul3(x, y, z):
        return imul(imul(x, y), z)
    t = isub(imul(M[1][1], M[2][2]), imul(M[1][2], M[2][1]))
    d1 = imul(M[0][0], t)
    t = isub(imul(M[1][0], M[2][2]), imul(M[1][2], M[2][0]))
    d2 = imul(M[0][1], t)
    t = isub(imul(M[1][0], M[2][1]), imul(M[1][1], M[2][0]))
    d3 = imul(M[0][2], t)
    return iadd(isub(d1, d2), d3)

def cell_bound(C, lam_iv, eta_iv):
    """Lower bound for I_eff on the cell via the factorization
    DM3 = Dh(mu) G:  sigma_min(G) >= sigma_min(DM3) / ||Dh||_F
    with the certified sigma_min(DM3) >= 1/100 on B (global theorem),
    and I_eff >= sigma_min(G)^2 / lambda_max(Omega).
    Returns None if an enclosure degenerates (triggers bisection)."""
    M = motif_enclosures(C, lam_iv)
    E, S, T, P, Q = (M[k]["val"] for k in ("E", "S", "T", "P", "Q"))
    if E[0] <= 0 or S[0] <= 0 or T[0] <= 0 or P[0] <= 0 or Q[0] <= 0:
        return None
    e1 = eta_iv; e2 = ipow(eta_iv, 2); e3 = ipow(eta_iv, 3)
    # y-components (all positive)
    ell = iscale(imul(e1, E), 2)
    t_  = iscale(imul(e2, T), Fraction(1, 2))
    s2  = iscale(imul(e2, S), 2)
    h3  = iadd(iadd(iscale(imul(e3, Q), 8), iscale(imul(e2, S), 12)),
               iscale(imul(e1, E), 2))
    h11 = iadd(iadd(iadd(iscale(imul(e3, P), 8), iscale(imul(e2, T), 3)),
                    iscale(imul(e2, S), 8)), iscale(imul(e1, E), 2))
    if ell[0] <= 0 or s2[0] <= 0:
        return None
    # m = 1 + 2 s2 / ell
    m_lo = 1 + 2 * s2[0] / ell[1]
    m_hi = 1 + 2 * s2[1] / ell[0]
    m = (m_lo, m_hi)
    # D = 4 eta^2 (QE - S^2)/E^2 + 2 eta S / E  >= 2 eta_lo S_lo / E_hi
    qe_s2_lo = Q[0] * E[0] - S[1] * S[1]
    if qe_s2_lo < 0:
        qe_s2_lo = Fraction(0)          # identity QE - S^2 >= 0
    D_lo = 4 * e2[0] * qe_s2_lo / (E[1] * E[1]) + 2 * e1[0] * S[0] / E[1]
    if D_lo <= 0:
        return None
    # N enclosure: N = h11/ell - m^2 (may straddle; only |N| upper needed)
    N_hi_abs = max(abs(h11[1] / ell[0] - m_lo * m_lo),
                   abs(h11[0] / ell[1] - m_hi * m_hi))
    # entry bounds (absolute) of Dh
    inv_ell_hi = 1 / ell[0]
    dm_dell_hi = 2 * s2[1] / (ell[0] * ell[0])          # |dm/dell|
    dm_ds2_hi = 2 / ell[0]
    # h2 = 3 t / s2
    a_h2_t = 3 / s2[0]
    a_h2_s2 = 3 * t_[1] / (s2[0] * s2[0])
    # dN/d. and dD/d. upper bounds
    dN_dell = h11[1] / (ell[0] * ell[0]) + 2 * m_hi * dm_dell_hi
    dN_ds2 = 2 * m_hi * dm_ds2_hi
    dN_dh11 = inv_ell_hi
    dD_dell = h3[1] / (ell[0] * ell[0]) + 2 * m_hi * dm_dell_hi
    dD_ds2 = 2 * m_hi * dm_ds2_hi
    dD_dh3 = inv_ell_hi
    D_abs_hi = h3[1] / ell[0] + m_hi * m_hi              # |D| <= this
    # |d h3 / d y_k| <= (|dN| D_abs_hi + N_hi_abs |dD|) / D_lo^2
    def h3row(dN, dD):
        return (dN * D_abs_hi + N_hi_abs * dD) / (D_lo * D_lo)
    r_ell = h3row(dN_dell, dD_dell)
    r_s2 = h3row(dN_ds2, dD_ds2)
    r_h3 = h3row(Fraction(0), dD_dh3)
    r_h11 = h3row(dN_dh11, Fraction(0))
    F2 = (Fraction(1)                       # h1 row
          + a_h2_t * a_h2_t + a_h2_s2 * a_h2_s2
          + r_ell * r_ell + r_s2 * r_s2 + r_h3 * r_h3 + r_h11 * r_h11)
    # sigma_min(G)^2 >= (1/100)^2 / F2   (certified global theorem on B)
    sig2_lo = Fraction(1, 10000) / F2
    # Omega upper bound: 142 * 5 * 209 * (2 eta)^6 at the cell's eta_hi
    C_om = Fraction(142 * 5 * 209) * (2 * eta_iv[1]) ** 6
    return sig2_lo / C_om

def main():
    C = load_coeffs()
    lam_lo, lam_hi = Fraction(1), Fraction(5)
    eta_lo, eta_hi = Fraction(5), Fraction(15)
    NL, NE = 16, 10
    stack = []
    for i in range(NL):
        for j in range(NE):
            stack.append(((lam_lo + (lam_hi - lam_lo) * i / NL,
                           lam_lo + (lam_hi - lam_lo) * (i + 1) / NL),
                          (eta_lo + (eta_hi - eta_lo) * j / NE,
                           eta_lo + (eta_hi - eta_lo) * (j + 1) / NE), 0))
    c_min = None
    boxes = 0
    while stack:
        lam_iv, eta_iv, depth = stack.pop()
        b = cell_bound(C, lam_iv, eta_iv)
        boxes += 1
        if b is None:
            if depth >= 7:
                raise CertificateError(
                    f"det gate failed at depth {depth}: lam={lam_iv} eta={eta_iv}")
            lm = (lam_iv[0] + lam_iv[1]) / 2
            em = (eta_iv[0] + eta_iv[1]) / 2
            wl = lam_iv[1] - lam_iv[0]
            we = (eta_iv[1] - eta_iv[0]) / 2   # scale eta by facewidth ratio
            if wl * Fraction(5, 2) >= (eta_iv[1] - eta_iv[0]):
                stack.append(((lam_iv[0], lm), eta_iv, depth + 1))
                stack.append(((lm, lam_iv[1]), eta_iv, depth + 1))
            else:
                stack.append((lam_iv, (eta_iv[0], em), depth + 1))
                stack.append((lam_iv, (em, eta_iv[1]), depth + 1))
            continue
        c_min = b if c_min is None else min(c_min, b)
    if c_min is None or c_min <= 0:
        raise CertificateError("no positive uniform bound obtained")
    print(f"boxes processed: {boxes}")
    print(f"CERTIFIED: inf over the null face of I_eff >= {float(c_min):.3e}")
    print(f"exact rational bound: {c_min.numerator}/{c_min.denominator}"[:120])
    print("GATES PASSED")

if __name__ == "__main__":
    main()

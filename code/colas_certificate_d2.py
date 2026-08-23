"""
colas_certificate_d2.py
=======================
Exact outward-rounded interval certificate for the TWO-DIMENSIONAL
hard-disc submodel.  The d=2 tangent-scale moment map equals the
d=1 map under the substitutions

    2*eta -> pi*eta^2 =: 2*eta_tilde     (so eta_tilde in [5,15] below)
    3/4   -> c2 = 1 - 3*sqrt(3)/(4*pi)   (the classical RGG clustering
                                          constant, enclosed by rational
                                          intervals of pi and sqrt(3))

so the entire coefficient-level machinery of colas_certificate_v2.py
applies verbatim with the geometric constant carried as an interval.
Certified conclusions on B2 = [1,5] x [5,15] x [0,3/20] in
(lambda, eta_tilde, psi):  global injectivity of the d=2 map
(psi=0 included) and uniform sigma_min(DM3) >= 1/100.
"""
from __future__ import annotations
import sys, time, math, os
from fractions import Fraction
from math import factorial

_BASE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, _BASE)
from colas_moments import build_coeffs, MOTIFS

PREC = 200
ONE = 1 << PREC
TRIM_TOL = Fraction(1, 1 << 96)     # coefficient-trim threshold -> error budget

# ---- d=2 geometric constant c2 = 1 - 3*sqrt(3)/(4*pi), enclosed ----------
from decimal import Decimal
def _dfrac(d):
    d = Decimal(d); return Fraction(int(d.scaleb(50)), 10**50)
_PI_LO = _dfrac("3.14159265358979323846264338327950288419716939937510")
_PI_HI = _PI_LO + Fraction(1, 10**49)
_S3_LO = _dfrac("1.73205080756887729352744634150587236694280525381038")
_S3_HI = _S3_LO + Fraction(1, 10**49)
C2_LO = 1 - 3*_S3_HI/(4*_PI_LO)
C2_HI = 1 - 3*_S3_LO/(4*_PI_HI)
assert Fraction(58,100) < C2_LO < C2_HI < Fraction(59,100)
CG = None  # set after IV is defined

# ----------------------------------------------------------------------
# fixed-point dyadic intervals
# ----------------------------------------------------------------------
def _fl(fr: Fraction) -> int:
    return (fr.numerator << PREC) // fr.denominator

def _ce(fr: Fraction) -> int:
    return -((-fr.numerator << PREC) // fr.denominator)

class IV:
    __slots__ = ("lo", "hi")
    def __init__(self, lo, hi): self.lo = lo; self.hi = hi
    @staticmethod
    def point(fr):
        fr = Fraction(fr); return IV(_fl(fr), _ce(fr))
    @staticmethod
    def box(lo, hi):
        return IV(_fl(Fraction(lo)), _ce(Fraction(hi)))
    def __add__(a, b): return IV(a.lo + b.lo, a.hi + b.hi)
    def __sub__(a, b): return IV(a.lo - b.hi, a.hi - b.lo)
    def __neg__(a):    return IV(-a.hi, -a.lo)
    def __mul__(a, b):
        p1 = a.lo*b.lo; p2 = a.lo*b.hi; p3 = a.hi*b.lo; p4 = a.hi*b.hi
        lo = min(p1,p2,p3,p4); hi = max(p1,p2,p3,p4)
        return IV(lo >> PREC, -((-hi) >> PREC))
    def sq(a):
        if a.lo >= 0: lo, hi = a.lo*a.lo, a.hi*a.hi
        elif a.hi <= 0: lo, hi = a.hi*a.hi, a.lo*a.lo
        else: lo, hi = 0, max(a.lo*a.lo, a.hi*a.hi)
        return IV(lo >> PREC, -((-hi) >> PREC))
    def __truediv__(a, b):
        if b.lo <= 0 <= b.hi: raise ZeroDivisionError
        cands = []
        for n in (a.lo, a.hi):
            for d in (b.lo, b.hi):
                num = n << PREC
                cands.append((num // d, -((-num) // d)))
        return IV(min(c[0] for c in cands), max(c[1] for c in cands))
    def widen_frac(a, e: Fraction):
        d = _ce(e); return IV(a.lo - d, a.hi + d)
    def scale_int(a, k):
        if k >= 0: return IV(a.lo*k, a.hi*k)
        return IV(a.hi*k, a.lo*k)
    def mag_frac(a):
        return Fraction(max(abs(a.lo), abs(a.hi)), ONE)
    @property
    def flo(self): return self.lo / ONE
    @property
    def fhi(self): return self.hi / ONE
    def frac_lo(self): return Fraction(self.lo, ONE)
    def frac_hi(self): return Fraction(self.hi, ONE)

IV0 = IV(0, 0)
IV1 = IV(ONE, ONE)
CG = IV(_fl(C2_LO), _ce(C2_HI))

def iv_sqrt_hi_frac(x_hi: Fraction) -> Fraction:
    """rational upper bound for sqrt(x_hi), x_hi >= 0."""
    n = (x_hi.numerator << (2 * PREC)) // x_hi.denominator
    s = math.isqrt(n) + 1
    return Fraction(s, ONE)


# ----------------------------------------------------------------------
# error-carrying polynomials in (dmu, psi): coeff[j] = list over dmu-powers
# value = sum_j psi^j sum_k c_{jk} dmu^k  +  err, |err| <= self.err
# ----------------------------------------------------------------------
class PE:
    __slots__ = ("c", "err", "sup")
    def __init__(self, c, err=Fraction(0), sup=None):
        self.c = c                    # dict j -> list[IV]
        self.err = err                # scalar error budget (Fraction, >=0)
        self.sup = sup                # crude sup bound on |value| over region

    @staticmethod
    def const(fr):
        return PE({0: [IV.point(fr)]}, Fraction(0), abs(Fraction(fr)))

    def copy(self):
        return PE({j: list(v) for j, v in self.c.items()}, self.err, self.sup)

    def __add__(a, b):
        c = {j: list(v) for j, v in a.c.items()}
        for j, v in b.c.items():
            if j in c:
                n = max(len(c[j]), len(v))
                arr = c[j] + [IV0]*(n-len(c[j]))
                for k in range(len(v)):
                    arr[k] = arr[k] + v[k]
                c[j] = arr
            else:
                c[j] = list(v)
        return PE(c, a.err + b.err, (a.sup or 0) + (b.sup or 0))

    def __sub__(a, b):
        c = {j: list(v) for j, v in a.c.items()}
        for j, v in b.c.items():
            if j in c:
                n = max(len(c[j]), len(v))
                arr = c[j] + [IV0]*(n-len(c[j]))
                for k in range(len(v)):
                    arr[k] = arr[k] - v[k]
                c[j] = arr
            else:
                c[j] = [-x for x in v]
        return PE(c, a.err + b.err, (a.sup or 0) + (b.sup or 0))

    def __neg__(a):
        return PE({j: [-x for x in v] for j, v in a.c.items()}, a.err, a.sup)

    def scale(a, fr):
        fr = Fraction(fr); t = IV.point(fr)
        return PE({j: [t * x for x in v] for j, v in a.c.items()},
                  abs(fr) * a.err, abs(fr) * (a.sup or 0))

    def scale_iv(a, t):
        mag = max(abs(t.frac_lo()), abs(t.frac_hi()))
        return PE({j: [t * x for x in v] for j, v in a.c.items()},
                  mag * a.err, mag * (a.sup or 0))

    def __mul__(a, b):
        c = {}
        for j1, v1 in a.c.items():
            for j2, v2 in b.c.items():
                j = j1 + j2
                n = len(v1) + len(v2) - 1
                arr = c.get(j)
                if arr is None:
                    arr = [IV0] * n; c[j] = arr
                elif len(arr) < n:
                    arr.extend([IV0] * (n - len(arr)))
                for k1, x1 in enumerate(v1):
                    if x1.lo == 0 and x1.hi == 0: continue
                    for k2, x2 in enumerate(v2):
                        arr[k1+k2] = arr[k1+k2] + x1 * x2
        sa = a.sup if a.sup is not None else Fraction(10**6)
        sb = b.sup if b.sup is not None else Fraction(10**6)
        err = sa * b.err + sb * a.err + a.err * b.err
        return PE(c, err, sa * sb)

    def trim(self, dmu_max: Fraction, psi_max: Fraction):
        """move tiny high-order coefficients into the error budget."""
        extra = Fraction(0)
        newc = {}
        for j, v in self.c.items():
            pj = psi_max ** j
            arr = []
            keep_to = -1
            mags = [x.mag_frac() for x in v]
            for k in range(len(v)):
                contrib = mags[k] * pj * (dmu_max ** k)
                if contrib > TRIM_TOL:
                    keep_to = k
            for k in range(keep_to + 1):
                arr.append(v[k])
            for k in range(keep_to + 1, len(v)):
                extra += mags[k] * pj * (dmu_max ** k)
            if arr:
                newc[j] = arr
        self.c = newc
        self.err += extra
        return self

    def eval(self, dmu: IV, psi: IV) -> IV:
        jmax = max(self.c.keys()) if self.c else 0
        acc = IV0
        for j in range(jmax, -1, -1):
            acc = acc * psi
            v = self.c.get(j)
            if v:
                h = v[-1]
                for x in reversed(v[:-1]):
                    h = h * dmu + x
                acc = acc + h
        return acc.widen_frac(self.err)


# ----------------------------------------------------------------------
# motif leaves per lambda cell
# ----------------------------------------------------------------------
N_TRUNC = 40
L_GLOBAL = Fraction(5)
EPS0 = Fraction(L_GLOBAL ** (N_TRUNC + 1), factorial(N_TRUNC + 1))
EPS1 = Fraction(L_GLOBAL ** N_TRUNC, factorial(N_TRUNC))
B_V = {2: 3, 3: 9, 4: 23}
SUP_VAL = Fraction(1)          # 0 <= motif <= 1
SUP_DL = Fraction(2)           # |d motif/d lambda| <= m * sup(uv e^{-l uv}) <= 3/e < 2
SUP_DP = {2: 3, 3: 9, 4: 23}   # |d motif/d psi| <= B_v

import pickle
_COEF_CACHE = os.path.join(_BASE, "coeffs_exact_N%d.pkl" % N_TRUNC)
if os.path.exists(_COEF_CACHE):
    with open(_COEF_CACHE, "rb") as _fh:
        _EXACT = pickle.load(_fh)
else:
    _EXACT = {k: build_coeffs(N_TRUNC, k) for k in "ESTPQ"}
    with open(_COEF_CACHE, "wb") as _fh:
        pickle.dump(_EXACT, _fh)

def shift_at_one(coeffs):
    c = list(coeffs); n = len(c)
    for i in range(n - 1):
        for j in range(n - 2, i - 1, -1):
            c[j] = c[j] + c[j + 1]
    return c

class Cell:
    """One lambda cell: holds motif leaves (val, dlam, dpsi) as PE polys
    in (dmu, psi), where lambda = lam0*(1+dmu)."""
    def __init__(self, lam_lo: Fraction, lam_hi: Fraction,
                 psi_hi: Fraction):
        self.lam_lo = lam_lo; self.lam_hi = lam_hi
        self.lam0 = (lam_lo + lam_hi) / 2
        self.dmu_max = (lam_hi - lam_lo) / 2 / self.lam0
        self.psi_hi = psi_hi
        self.leaf = {}
        inv0 = Fraction(1) / self.lam0
        for k in "ESTPQ":
            a0, a1, a2 = _EXACT[k]
            v = MOTIFS[k]["v"]; m = len(MOTIFS[k]["edges"])
            c1 = (1 + EPS0) ** (m - 1)
            c2 = (1 + EPS0) ** max(m - 2, 0)
            err_v = (2 ** v) * m * EPS0 * c1
            err_l = (2 ** v) * m * (EPS1 * c1 + (m - 1) * EPS0 * c2)
            err_p = B_V[v] * m * EPS0 * c1
            rows_v, rows_l, rows_p = {}, {}, {}
            for j, arr in enumerate((a0, a1, a2)):
                if all(x == 0 for x in arr):
                    continue
                pw = Fraction(1)
                resc = []
                for cc in arr:
                    resc.append(IV.point(cc * pw)); pw *= self.lam0
                s = shift_at_one(resc)
                rows_v[j] = s
                d = [s[t].scale_int(t) for t in range(1, len(s))]
                iv0 = IV.point(inv0)
                rows_l[j] = [iv0 * x for x in d]
                if j >= 1:
                    prev = rows_p.get(j - 1, [IV0] * len(s))
                    n = max(len(prev), len(s))
                    prev = prev + [IV0] * (n - len(prev))
                    for t in range(len(s)):
                        prev[t] = prev[t] + s[t].scale_int(j)
                    rows_p[j - 1] = prev
            val = PE(rows_v, err_v, SUP_VAL).trim(self.dmu_max, psi_hi)
            dl = PE(rows_l, err_l, SUP_DL).trim(self.dmu_max, psi_hi)
            dp = PE(rows_p, err_p, SUP_DP[v]).trim(self.dmu_max, psi_hi)
            self.leaf[k] = (val, dl, dp)

    # ---------------- G1 composite polynomials ------------------------
    def build_g1(self):
        E, El, Ep = self.leaf["E"]; S, Sl, Sp = self.leaf["S"]
        T, Tl, Tp = self.leaf["T"]; P, Pl, Pp = self.leaf["P"]
        Q, Ql, Qp = self.leaf["Q"]
        tr = lambda p: p.trim(self.dmu_max, self.psi_hi)
        A  = tr(P*E - S*S);  B = tr(Q*E - S*S)
        Al = tr(Pl*E + P*El - (S*Sl).scale(2))
        Ap = tr(Pp*E + P*Ep - (S*Sp).scale(2))
        Bl = tr(Ql*E + Q*El - (S*Sl).scale(2))
        Bp = tr(Qp*E + Q*Ep - (S*Sp).scale(2))
        E2 = tr(E*E)
        TE2 = tr(T*E2); SE2 = tr(S*E2)
        TE2l = tr(Tl*E2 + (T*E*El).scale(2))
        TE2p = tr(Tp*E2 + (T*E*Ep).scale(2))
        SE2l = tr(Sl*E2 + (S*E*El).scale(2))
        SE2p = tr(Sp*E2 + (S*E*Ep).scale(2))
        NC  = tr(Tl*S - T*Sl)          # C_lambda * (4/3) * S^2
        NCp = tr(Tp*S - T*Sp)          # C_psi    * (4/3) * S^2
        c34 = None  # geometric constant CG used instead
        # num = ell*A + (3/4)TE2 ; den = ell*B + SE2
        # N1 = num_l*den - num*den_l  (quadratic in ell)
        N1 = [tr(TE2l.scale_iv(CG)*SE2 - TE2.scale_iv(CG)*SE2l),
              tr(Al*SE2 + TE2l.scale_iv(CG)*B - A*SE2l - TE2.scale_iv(CG)*Bl),
              tr(Al*B - A*Bl)]
        N2 = [tr(TE2p.scale_iv(CG)*SE2 - TE2.scale_iv(CG)*SE2p),
              tr(Ap*SE2 + TE2p.scale_iv(CG)*B - A*SE2p - TE2.scale_iv(CG)*Bp),
              tr(Ap*B - A*Bp)]
        # DeltaN = N2*NC - N1*NCp  (quadratic in ell); sign of Delta
        DN = [tr(N2[i]*NC - N1[i]*NCp) for i in range(3)]
        den_c = [SE2, B]               # den' = ell*B + SE2
        num_c = [TE2.scale_iv(CG), A]
        return dict(S=S, NC=NC, NCp=NCp, DN=DN, den=den_c, num=num_c,
                    S2=tr(S*S))

    # ---------------- G2 composite polynomials ------------------------
    def build_g2(self):
        E, El, Ep = self.leaf["E"]; S, Sl, Sp = self.leaf["S"]
        T, Tl, Tp = self.leaf["T"]; P, Pl, Pp = self.leaf["P"]
        Q, Ql, Qp = self.leaf["Q"]
        tr = lambda p: p.trim(self.dmu_max, self.psi_hi)
        c34 = None  # geometric constant CG used instead
        A  = tr(P*E - S*S);  B = tr(Q*E - S*S)
        Al = tr(Pl*E + P*El - (S*Sl).scale(2))
        Ap = tr(Pp*E + P*Ep - (S*Sp).scale(2))
        Bl = tr(Ql*E + Q*El - (S*Sl).scale(2))
        Bp = tr(Qp*E + Q*Ep - (S*Sp).scale(2))
        TE = tr(T*E); SE = tr(S*E)
        TEl = tr(Tl*E + T*El); TEp = tr(Tp*E + T*Ep)
        SEl = tr(Sl*E + S*El); SEp = tr(Sp*E + S*Ep)
        NC  = tr(Tl*S - T*Sl)
        NCp = tr(Tp*S - T*Sp)
        # num = 2eta*A + (3/4)TE ; den = 2eta*B + SE   (eta-linear)
        # r_x = Nr_x / den^2 with (eta-quadratic numerators):
        # Nr_l = num_l*den - num*den_l, similarly _p; Nr_e = num_e*den - num*den_e
        # with num_e = 2A, den_e = 2B  (d/d eta of 2eta = 2)
        num = [TE.scale_iv(CG), A.scale(2)]      # [eta^0, eta^1]
        den = [SE, B.scale(2)]
        num_l = [TEl.scale_iv(CG), Al.scale(2)]
        num_p = [TEp.scale_iv(CG), Ap.scale(2)]
        den_l = [SEl, Bl.scale(2)]
        den_p = [SEp, Bp.scale(2)]
        num_e = [A.scale(2)]                   # eta^0
        den_e = [B.scale(2)]
        def polymul_eta(u, v):
            out = [None] * (len(u) + len(v) - 1)
            for i, ui in enumerate(u):
                for j, vj in enumerate(v):
                    t = ui * vj
                    out[i+j] = t if out[i+j] is None else tr(out[i+j] + t)
            return [tr(x) for x in out]
        def polysub_eta(u, v):
            n = max(len(u), len(v))
            zero = PE({}, Fraction(0), Fraction(0))
            u = u + [zero] * (n - len(u)); v = v + [zero] * (n - len(v))
            return [tr(u[i] - v[i]) for i in range(n)]
        Nr_l = polysub_eta(polymul_eta(num_l, den), polymul_eta(num, den_l))
        Nr_p = polysub_eta(polymul_eta(num_p, den), polymul_eta(num, den_p))
        Nr_e = polysub_eta(polymul_eta(num_e, den), polymul_eta(num, den_e))
        # ell-row entries (polynomials in eta):
        ell_l = [PE({}, Fraction(0), Fraction(0)), El.scale(2)]
        ell_e = [E.scale(2)]
        ell_p = [PE({}, Fraction(0), Fraction(0)), Ep.scale(2)]
        # det numerator: det = [ -NC*(ell_e*Nr_p - ell_p*Nr_e)
        #                        -NCp*(ell_l*Nr_e - ell_e*Nr_l) ] * (3/4) / (S^2 den^2)
        m_ep = polysub_eta(polymul_eta(ell_e, Nr_p), polymul_eta(ell_p, Nr_e))
        m_le = polysub_eta(polymul_eta(ell_l, Nr_e), polymul_eta(ell_e, Nr_l))
        DNdet = polysub_eta([x.scale(-1) for x in polymul_eta([NC], m_ep)],
                            polymul_eta([NCp], m_le))
        DNdet = [x.scale_iv(CG) for x in DNdet]
        # 2x2 minor numerators for Lambda^2 Frobenius:
        # rows: R1=(ell_l, ell_e, ell_p) [poly], R2=(3/4)(NC,0,NCp)/S^2,
        #       R3=(Nr_l,Nr_e,Nr_p)/den^2
        pairs = [(0,1),(0,2),(1,2)]
        R1 = [ell_l, ell_e, ell_p]
        R2n = [[NC.scale_iv(CG)], [PE({}, Fraction(0), Fraction(0))], [NCp.scale_iv(CG)]]
        R3n = [Nr_l, Nr_e, Nr_p]
        M12, M13, M23 = [], [], []
        for (c1, c2) in pairs:
            M12.append(polysub_eta(polymul_eta(R1[c1], R2n[c2]),
                                   polymul_eta(R1[c2], R2n[c1])))
            M13.append(polysub_eta(polymul_eta(R1[c1], R3n[c2]),
                                   polymul_eta(R1[c2], R3n[c1])))
            M23.append(polysub_eta(polymul_eta(R2n[c1], R3n[c2]),
                                   polymul_eta(R2n[c2], R3n[c1])))
        return dict(S=S, S2=tr(S*S), den=den, DNdet=DNdet,
                    M12=M12, M13=M13, M23=M23,
                    R1=R1, R2n=R2n, R3n=R3n, NC=NC, NCp=NCp)


def eval_eta(poly_eta, dmu: IV, psi: IV, eta: IV) -> IV:
    acc = IV0
    for coef in reversed(poly_eta):
        acc = acc * eta + coef.eval(dmu, psi)
    return acc


# ----------------------------------------------------------------------
# gate constants and drivers
# ----------------------------------------------------------------------
G_DET = Fraction(-1, 64)
G_SMIN = Fraction(1, 100)
LAM_LO, LAM_HI = Fraction(1), Fraction(5)
PSI_LO, PSI_HI = Fraction(0), Fraction(3, 20)
ETA_LO, ETA_HI = Fraction(5), Fraction(15)

N_CELLS = 32          # fixed lambda cells of width 1/8 over [1,5]
CELL_W = (LAM_HI - LAM_LO) / N_CELLS

def cell_bounds(idx):
    return LAM_LO + CELL_W * idx, LAM_LO + CELL_W * (idx + 1)

class CellCache:
    """Composite polynomials built ONCE per fixed lambda cell; boxes with
    lambda-bounds inside a cell reuse the cell's polynomials and only
    shrink the dmu evaluation interval."""
    def __init__(self, psi_hi, builder):
        self.psi_hi = psi_hi; self.builder = builder; self.d = {}
    def get(self, idx):
        if idx not in self.d:
            lo, hi = cell_bounds(idx)
            cell = Cell(lo, hi, self.psi_hi)
            self.d[idx] = (getattr(cell, self.builder)(), cell)
        return self.d[idx]


def certify_E_range(cache, max_boxes=400000):
    stack = []
    for i in range(N_CELLS):
        lo, hi = cell_bounds(i)
        stack.append((i, lo, hi, PSI_LO, PSI_HI, 0))
    elo = ehi = None
    processed = 0
    TOL = Fraction(1, 500)
    while stack:
        idx, lam_lo, lam_hi, ps_lo, ps_hi, depth = stack.pop()
        processed += 1
        if processed > max_boxes: raise RuntimeError("E budget")
        _, cell = cache.get(idx)
        Ep = cell.leaf["E"][0]
        dmu = IV.box((lam_lo - cell.lam0) / cell.lam0, (lam_hi - cell.lam0) / cell.lam0)
        val = Ep.eval(dmu, IV.box(ps_lo, ps_hi))
        wid = val.frac_hi() - val.frac_lo()
        if wid > TOL and depth < 18:
            if (lam_hi - lam_lo) / CELL_W >= (ps_hi - ps_lo) / PSI_HI:
                mid = (lam_lo + lam_hi) / 2
                stack += [(idx, lam_lo, mid, ps_lo, ps_hi, depth + 1),
                          (idx, mid, lam_hi, ps_lo, ps_hi, depth + 1)]
            else:
                mid = (ps_lo + ps_hi) / 2
                stack += [(idx, lam_lo, lam_hi, ps_lo, mid, depth + 1),
                          (idx, lam_lo, lam_hi, mid, ps_hi, depth + 1)]
            continue
        lo, hi = val.frac_lo(), val.frac_hi()
        elo = lo if elo is None else min(elo, lo)
        ehi = hi if ehi is None else max(ehi, hi)
    return elo, ehi, processed


def certify_G1(cache, ell_lo, ell_hi, max_boxes=800000, progress=True):
    stack = []
    for i in range(N_CELLS):
        lo, hi = cell_bounds(i)
        stack.append((i, lo, hi, PSI_LO, PSI_HI, ell_lo, ell_hi, 0))
    processed = leaves = 0
    st = dict(cl_lo=None, cp_lo=None, denp_lo=None, slope_lo=None,
              slope_hi=None, delta_lo=None)
    t0 = time.time()
    while stack:
        idx, lam_lo, lam_hi, ps_lo, ps_hi, el_lo, el_hi, depth = stack.pop()
        processed += 1
        if processed > max_boxes: raise RuntimeError("G1 budget")
        if progress and processed % 5000 == 0:
            print("  [G1] boxes=%d stack=%d ok=%d t=%.0fs"
                  % (processed, len(stack), leaves, time.time() - t0),
                  file=sys.stderr, flush=True)
        g, cell = cache.get(idx)
        dmu = IV.box((lam_lo - cell.lam0) / cell.lam0,
                     (lam_hi - cell.lam0) / cell.lam0)
        psi = IV.box(ps_lo, ps_hi); ell = IV.box(el_lo, el_hi)
        Sv = g["S"].eval(dmu, psi)
        NCv = g["NC"].eval(dmu, psi)
        NCpv = g["NCp"].eval(dmu, psi)
        denv = eval_eta(g["den"], dmu, psi, ell)     # den' (ell plays eta role)
        DNv = eval_eta(g["DN"], dmu, psi, ell)
        ok = (Sv.lo > 0 and denv.lo > 0 and NCv.lo > 0 and NCpv.lo > 0
              and DNv.lo > 0)
        if ok:
            leaves += 1
            S2v = g["S2"].eval(dmu, psi)
            cl = (CG * NCv) / S2v
            cp = (CG * NCpv) / S2v
            dlt = (CG * DNv) / (denv.sq() * S2v)
            slope = dlt / cl
            for key, val in (("cl_lo", cl.frac_lo()), ("cp_lo", cp.frac_lo()),
                             ("denp_lo", denv.frac_lo()), ("delta_lo", dlt.frac_lo())):
                st[key] = val if st[key] is None else min(st[key], val)
            st["slope_lo"] = slope.frac_lo() if st["slope_lo"] is None else min(st["slope_lo"], slope.frac_lo())
            st["slope_hi"] = slope.frac_hi() if st["slope_hi"] is None else max(st["slope_hi"], slope.frac_hi())
            continue
        if depth >= 30:
            raise RuntimeError("G1 depth cap at box (%s,%s)x(%s,%s)x(%s,%s)"
                               % (lam_lo, lam_hi, ps_lo, ps_hi, el_lo, el_hi))
        wl = (lam_hi - lam_lo) / CELL_W
        wp = (ps_hi - ps_lo) / (PSI_HI - PSI_LO)
        we = (el_hi - el_lo) / (ell_hi - ell_lo)
        wmax = max(wl, wp, we)
        if wmax == wl:
            mid = (lam_lo + lam_hi) / 2
            stack += [(idx, lam_lo, mid, ps_lo, ps_hi, el_lo, el_hi, depth + 1),
                      (idx, mid, lam_hi, ps_lo, ps_hi, el_lo, el_hi, depth + 1)]
        elif wmax == wp:
            mid = (ps_lo + ps_hi) / 2
            stack += [(idx, lam_lo, lam_hi, ps_lo, mid, el_lo, el_hi, depth + 1),
                      (idx, lam_lo, lam_hi, mid, ps_hi, el_lo, el_hi, depth + 1)]
        else:
            mid = (el_lo + el_hi) / 2
            stack += [(idx, lam_lo, lam_hi, ps_lo, ps_hi, el_lo, mid, depth + 1),
                      (idx, lam_lo, lam_hi, ps_lo, ps_hi, mid, el_hi, depth + 1)]
    return processed, leaves, st


def certify_G2(cache, max_boxes=1200000, progress=True):
    stack = []
    for i in range(N_CELLS):
        lo, hi = cell_bounds(i)
        stack.append((i, lo, hi, ETA_LO, ETA_HI, PSI_LO, PSI_HI, 0))
    processed = leaves = 0
    st = dict(det_lo=None, det_hi=None, frob2_hi=None, smin_lo=None)
    t0 = time.time()
    while stack:
        idx, lam_lo, lam_hi, et_lo, et_hi, ps_lo, ps_hi, depth = stack.pop()
        processed += 1
        if processed > max_boxes: raise RuntimeError("G2 budget")
        if progress and processed % 5000 == 0:
            print("  [G2] boxes=%d stack=%d ok=%d t=%.0fs"
                  % (processed, len(stack), leaves, time.time() - t0),
                  file=sys.stderr, flush=True)
        g, cell = cache.get(idx)
        dmu = IV.box((lam_lo - cell.lam0) / cell.lam0,
                     (lam_hi - cell.lam0) / cell.lam0)
        psi = IV.box(ps_lo, ps_hi); eta = IV.box(et_lo, et_hi)
        Sv = g["S"].eval(dmu, psi)
        S2v = g["S2"].eval(dmu, psi)
        denv = eval_eta(g["den"], dmu, psi, eta)
        ok = False
        if Sv.lo > 0 and denv.lo > 0 and S2v.lo > 0:
            DNv = eval_eta(g["DNdet"], dmu, psi, eta)
            det = DNv / (S2v * denv.sq())
            if det.frac_hi() < G_DET:
                # Lambda^2 Frobenius^2
                f2 = IV0
                den2 = denv.sq(); S4 = S2v.sq(); d4 = den2.sq()
                for Mv in g["M12"]:
                    f2 = f2 + eval_eta(Mv, dmu, psi, eta).sq() / S4
                for Mv in g["M13"]:
                    f2 = f2 + eval_eta(Mv, dmu, psi, eta).sq() / d4
                for Mv in g["M23"]:
                    f2 = f2 + eval_eta(Mv, dmu, psi, eta).sq() / (S4 * d4)
                smin_lo = (-det.frac_hi()) / iv_sqrt_hi_frac(f2.frac_hi())
                if smin_lo > G_SMIN:
                    ok = True
        if ok:
            leaves += 1
            # plain Frobenius^2 for reporting
            fr2 = IV0
            for e in g["R1"]:
                fr2 = fr2 + eval_eta(e, dmu, psi, eta).sq()
            for e in g["R2n"]:
                fr2 = fr2 + (eval_eta(e, dmu, psi, eta) / S2v).sq()
            for e in g["R3n"]:
                fr2 = fr2 + (eval_eta(e, dmu, psi, eta) / denv.sq()).sq()
            st["det_lo"] = det.frac_lo() if st["det_lo"] is None else min(st["det_lo"], det.frac_lo())
            st["det_hi"] = det.frac_hi() if st["det_hi"] is None else max(st["det_hi"], det.frac_hi())
            st["frob2_hi"] = fr2.frac_hi() if st["frob2_hi"] is None else max(st["frob2_hi"], fr2.frac_hi())
            st["smin_lo"] = smin_lo if st["smin_lo"] is None else min(st["smin_lo"], smin_lo)
            continue
        if depth >= 30:
            raise RuntimeError("G2 depth cap at (%s,%s)x(%s,%s)x(%s,%s)"
                               % (lam_lo, lam_hi, et_lo, et_hi, ps_lo, ps_hi))
        wl = (lam_hi - lam_lo) / CELL_W
        we = (et_hi - et_lo) / (ETA_HI - ETA_LO)
        wp = (ps_hi - ps_lo) / (PSI_HI - PSI_LO)
        wmax = max(wl, we, wp)
        if wmax == wl:
            mid = (lam_lo + lam_hi) / 2
            stack += [(idx, lam_lo, mid, et_lo, et_hi, ps_lo, ps_hi, depth + 1),
                      (idx, mid, lam_hi, et_lo, et_hi, ps_lo, ps_hi, depth + 1)]
        elif wmax == we:
            mid = (et_lo + et_hi) / 2
            stack += [(idx, lam_lo, lam_hi, et_lo, mid, ps_lo, ps_hi, depth + 1),
                      (idx, lam_lo, lam_hi, mid, et_hi, ps_lo, ps_hi, depth + 1)]
        else:
            mid = (ps_lo + ps_hi) / 2
            stack += [(idx, lam_lo, lam_hi, et_lo, et_hi, ps_lo, mid, depth + 1),
                      (idx, lam_lo, lam_hi, et_lo, et_hi, mid, ps_hi, depth + 1)]
    return processed, leaves, st


def selftest():
    a = IV.box(Fraction(1, 3), Fraction(1, 2)); b = IV.box(Fraction(-2), Fraction(3))
    c = a * b
    assert c.flo <= -1.0 and c.fhi >= 1.5
    d = IV.box(2, 3) / IV.box(4, 5)
    assert d.flo <= 0.4 and d.fhi >= 0.75
    # speck-box reproduction of v1 (single cell, one box; poly-level composites)
    lam_lo, lam_hi = Fraction(299, 100), Fraction(301, 100)
    cell = Cell(lam_lo, lam_hi, Fraction(87, 1000))
    g2 = cell.build_g2()
    dmu = IV.box((lam_lo - cell.lam0) / cell.lam0, (lam_hi - cell.lam0) / cell.lam0)
    psi = IV.box(Fraction(2, 25), Fraction(87, 1000))
    eta = IV.box(Fraction(89, 10), Fraction(91, 10))
    S2v = g2["S2"].eval(dmu, psi)
    denv = eval_eta(g2["den"], dmu, psi, eta)
    DNv = eval_eta(g2["DNdet"], dmu, psi, eta)
    det = DNv / (S2v * denv.sq())
    print("selftest speck box (single box, no subdivision):")
    print("  det J in [%.6f, %.6f]   (v1, 1100 boxes: [-0.078436, -0.055595])"
          % (det.flo, det.fhi))
    f2 = IV0
    den2 = denv.sq(); S4 = S2v.sq(); d4 = den2.sq()
    for Mv in g2["M12"]: f2 = f2 + eval_eta(Mv, dmu, psi, eta).sq() / S4
    for Mv in g2["M13"]: f2 = f2 + eval_eta(Mv, dmu, psi, eta).sq() / d4
    for Mv in g2["M23"]: f2 = f2 + eval_eta(Mv, dmu, psi, eta).sq() / (S4 * d4)
    print("  sigma_min >= %.6f       (v1 bound: 0.017657)"
          % (-det.fhi / float(iv_sqrt_hi_frac(f2.frac_hi()))))
    g1 = cell.build_g1()
    ell = IV.box(Fraction(78, 10), Fraction(81, 10))
    NCv = g1["NC"].eval(dmu, psi); S2 = g1["S2"].eval(dmu, psi)
    cl = NCv.scale_int(3) / S2.scale_int(4)
    denp = eval_eta(g1["den"], dmu, psi, ell)
    DN = eval_eta(g1["DN"], dmu, psi, ell)
    dlt = DN.scale_int(3) / (denp.sq() * S2).scale_int(4)
    print("  C_lambda in [%.6f, %.6f] (v1 min: 0.061628; float center 0.062318)"
          % (cl.flo, cl.fhi))
    print("  slope in [%.6f, %.6f]    (float center 1.211562)"
          % ((dlt / cl).flo, (dlt / cl).fhi))
    print("selftest done")


def main():
    if "--selftest" in sys.argv:
        t0 = time.time(); selftest(); print("t=%.1fs" % (time.time() - t0)); return
    t0 = time.time()
    print("Certificate d2: two-dimensional hard-disc submodel")
    print("B2 = [1,5] x [5,15] x [0,3/20]  (lambda, eta_tilde = pi eta^2/2, psi);  PREC=%d, N=%d" % (PREC, N_TRUNC))
    print("geometric constant c2 = 1 - 3 sqrt(3)/(4 pi) enclosed in [%.12f, %.12f]" % (CG.flo, CG.fhi))
    cache_g1 = CellCache(PSI_HI, "build_g1")
    print("--- Stage 0: E-range on [1,5]x[0,3/20] ---", flush=True)
    elo, ehi, nb = certify_E_range(cache_g1)
    print("E in [%.9f, %.9f]   (%d boxes)" % (float(elo), float(ehi), nb), flush=True)
    ell_lo, ell_hi = 10 * elo, 30 * ehi
    print("ell-image of B inside [%.6f, %.6f]" % (float(ell_lo), float(ell_hi)))
    print("--- Stage 1 (G1): injectivity gates ---", flush=True)
    p1, l1, s1 = certify_G1(cache_g1, ell_lo, ell_hi)
    # explicit theorem-level gates (active under python -O as well)
    if not (elo >= Fraction(2034, 10000) and ehi <= Fraction(56805, 100000)):
        raise RuntimeError("E-range gate failed")
    if not s1["cl_lo"] > Fraction(1, 64):
        raise RuntimeError("gate C_lambda > 1/64 failed")
    if not s1["cp_lo"] > Fraction(1, 16):
        raise RuntimeError("gate C_psi > 1/16 failed")
    if not s1["delta_lo"] > Fraction(1, 800):
        raise RuntimeError("gate Delta > 1/800 failed")
    if not s1["denp_lo"] > Fraction(1, 512):
        raise RuntimeError("gate den' > 1/512 failed")
    print("G1 PASS with %d boxes (%d leaves)" % (p1, l1))
    print("  C_lambda >= %.6f  (gate 1/64) ; C_psi >= %.6f  (gate 1/16)"
          % (float(s1["cl_lo"]), float(s1["cp_lo"])))
    print("  Delta    >= %.6f  (gate 1/800); den' >= %.6f  (gate 1/512)"
          % (float(s1["delta_lo"]), float(s1["denp_lo"])))
    print("  compensated slope in [%.6f, %.6f]" % (float(s1["slope_lo"]), float(s1["slope_hi"])))
    del cache_g1
    cache_g2 = CellCache(PSI_HI, "build_g2")
    print("--- Stage 2 (G2): uniform rank on B ---", flush=True)
    p2, l2, s2 = certify_G2(cache_g2)
    if not s2["det_hi"] < Fraction(-1, 64):
        raise RuntimeError("gate det < -1/64 failed")
    if not s2["smin_lo"] > Fraction(1, 100):
        raise RuntimeError("gate sigma_min > 1/100 failed")
    print("G2 PASS with %d boxes (%d leaves)" % (p2, l2))
    print("  det DM3 in [%.6f, %.6f]  (< -1/64 uniformly)" % (float(s2["det_lo"]), float(s2["det_hi"])))
    print("  sup ||DM3||_F^2 <= %.5f" % float(s2["frob2_hi"]))
    print("  certified sigma_min lower bound = %.6f  (theorem constant 1/100)" % float(s2["smin_lo"]))
    print("PASS: d=2 hard-disc M3 injective on B2 (psi=0 included); sigma_min >= 1/100 on B2")
    print("total time %.1fs" % (time.time() - t0), file=sys.stderr)


if __name__ == "__main__":
    main()

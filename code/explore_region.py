"""Map the moment map over the enlarged region to choose certifiable boxes."""
import numpy as np, pickle, sys
import os as _os0
sys.path.insert(0, _os0.path.dirname(_os0.path.abspath(__file__)))
import os as _os
_BASE = _os.environ.get("COLAS_BASE") or _os.path.dirname(
    _os.path.dirname(_os.path.abspath(__file__)))
def _p(*parts):
    q = _os.path.join(_BASE, *parts)
    _os.makedirs(_os.path.dirname(q), exist_ok=True) if parts[:1] in (("runs",), ("figures",), ("revision",)) else None
    return q

from colas_moments import MomentMap

mm = MomentMap(N=60, cache=_p("runs", "coeffs_N60.pkl"))

# ---------------- (lambda, psi) grid: C_lambda, motifs, C, slope over ell ------
lam = np.linspace(0.4, 6.0, 141)
psi = np.linspace(0.0, 0.16, 81)
L, PS = np.meshgrid(lam, psi, indexing="ij")

red_any = {}
ells = np.array([1.0, 2.0, 4.0, 8.0, 12.0, 16.0, 20.0])
slope_min = np.full(L.shape, np.inf)
slope_max = np.full(L.shape, -np.inf)
Delta_min = np.full(L.shape, np.inf)
for ell in ells:
    red = mm.reduced(L, PS, np.full(L.shape, ell))
    slope_min = np.minimum(slope_min, red["slope"])
    slope_max = np.maximum(slope_max, red["slope"])
    Delta_min = np.minimum(Delta_min, red["Delta"])
red0 = mm.reduced(L, PS, np.full(L.shape, 8.0))
C_l = red0["C_l"]

# targets on the same grid at a reference eta (eta only scales ell / shifts r)
t = mm.targets(L, np.full(L.shape, 9.0), PS)
C = t["C"]; E = t["motifs"]["E"]

# ---------------- sigma_min over (lambda, eta, psi) ----------------------------
lam3 = np.linspace(0.5, 6.0, 56)
eta3 = np.linspace(4.0, 16.0, 25)
psi3 = np.linspace(0.0, 0.16, 33)
L3, H3, P3 = np.meshgrid(lam3, eta3, psi3, indexing="ij")
t3 = mm.targets(L3, H3, P3)
J3 = t3["J"]
sv = np.linalg.svd(J3, compute_uv=False)
smin = sv[..., -1]
detJ3 = np.linalg.det(J3)

out = dict(lam=lam, psi=psi, C_l=C_l, slope_min=slope_min, slope_max=slope_max,
           Delta_min=Delta_min, C=C, E=E,
           lam3=lam3, eta3=eta3, psi3=psi3, smin=smin, detJ3=detJ3,
           r_grid=t["r"], ell_at_eta9=t["ell"])
with open(_p("runs", "explore.pkl"), "wb") as fh:
    pickle.dump(out, fh)

print("C_lambda: min over grid = %.6f at lam=%.2f psi=%.3f" %
      (C_l.min(), L.flat[C_l.argmin()], PS.flat[C_l.argmin()]))
print("Delta:    min over grid x ells = %.6f at lam=%.2f psi=%.3f" %
      (Delta_min.min(), L.flat[Delta_min.argmin()], PS.flat[Delta_min.argmin()]))
print("slope:    range [%.4f, %.4f]" % (slope_min.min(), slope_max.max()))
print("sigma_min over (lam,eta,psi) box [0.5,6]x[4,16]x[0,0.16]:")
print("  min = %.6f at lam=%.2f eta=%.2f psi=%.3f" %
      (smin.min(), L3.flat[smin.argmin()], H3.flat[smin.argmin()], P3.flat[smin.argmin()]))
for lo, hi in [(1.0, 5.0), (1.5, 5.0), (2.0, 5.0)]:
    m = (L3 >= lo) & (L3 <= hi) & (H3 >= 5) & (H3 <= 15)
    print("  min over lam in [%.1f,%.1f], eta in [5,15], psi in [0,0.16]: %.6f" %
          (lo, hi, smin[m].min()))
print("detJ sign: max over grid = %.6f (should be <0 everywhere if no sign change)" % detJ3.max())
print("C range over (lam,psi) grid: [%.4f, %.4f]" % (C.min(), C.max()))
print("E range: [%.4f, %.4f] -> ell=2*eta*E over eta in [5,15]: [%.2f, %.2f]" %
      (E.min(), E.max(), 2 * 5 * E.min(), 2 * 15 * E.max()))
print("r range at eta=9: [%.4f, %.4f]" % (t["r"].min(), t["r"].max()))

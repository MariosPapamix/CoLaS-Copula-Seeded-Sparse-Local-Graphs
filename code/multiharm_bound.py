"""Worst-case effect of the multi-harmonic e2-perturbation on the
one-harmonic inversion (SI S5, Corollary): scales the psi^2 components of
the four-vertex motifs by kappa = 1 + 2 e2/psi^2 <= 2 and re-inverts."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from colas_moments import MomentMap
from colas_sim import Estimator
mm = MomentMap(N=60); est = Estimator(mm)
def perturbed(lam, eta, psi, kappa=2.0):
    M = {}
    for k in "ESTPQ":
        a0, a1, a2 = mm.fc[k]
        s = np.arange(len(a0)); L = lam ** s
        kap = kappa if k in "PQ" else 1.0
        M[k] = float((L * (a0 + psi*a1 + kap*psi*psi*a2)).sum())
    E,S,T,P,Q = (M[k] for k in "ESTPQ")
    return np.array([2*eta*E, 0.75*T/S,
                     (2*eta*(P*E-S*S)+0.75*T*E)/(2*eta*(Q*E-S*S)+S*E)])
worst = 0.0
for th in [(1.5,6,0.05),(3,9,0.0835),(4.5,14,0.14),(2,12,0.10),(3,9,0.15)]:
    th_hat, ok = est.solve_mom(perturbed(*th), x0=th)
    rel = abs(th_hat[2]-th[2])/th[2]
    worst = max(worst, rel)
    print("theta=%s -> psi_hat=%.5f (rel err %.2f%%)" % (th, th_hat[2], 100*rel))
print("worst relative psi_tot error: %.2f%%" % (100*worst))
assert worst < 0.03

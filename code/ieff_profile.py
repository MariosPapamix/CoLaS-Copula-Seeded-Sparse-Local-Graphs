"""Numerical profile of the efficient boundary information I_eff along
the null face (honest magnitudes to accompany the conservative
certificate).  I_eff = 1 / [(G^T Omega^-1 G)^{-1}]_{psi,psi} with G
from the exact coefficient tables and Omega estimated by simulation of
the five-moment vector at n = 10^4.  Writes runs/ieff_profile.json.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from colas_moments import MomentMap
import colas_sim as cs

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
N_SIM, REPS = 10000, 160

def G_at(mm, lam, eta):
    M = mm.all_motifs(np.array(lam), np.array(0.0))
    E, El, Ep = M["E"]; S, Sl, Sp = M["S"]; T, Tl, Tp = M["T"]
    P, Pl, Pp = M["P"]; Q, Ql, Qp = M["Q"]
    return np.array([
        [2*eta*El, 2*E, 2*eta*Ep],
        [0.5*eta**2*Tl, eta*T, 0.5*eta**2*Tp],
        [2*eta**2*Sl, 4*eta*S, 2*eta**2*Sp],
        [8*eta**3*Ql+12*eta**2*Sl+2*eta*El, 24*eta**2*Q+24*eta*S+2*E,
         8*eta**3*Qp+12*eta**2*Sp+2*eta*Ep],
        [8*eta**3*Pl+3*eta**2*Tl+8*eta**2*Sl+2*eta*El,
         24*eta**2*P+6*eta*T+16*eta*S+2*E,
         8*eta**3*Pp+3*eta**2*Tp+8*eta**2*Sp+2*eta*Ep]], float)

def omega_hat(lam, eta, rng):
    Y = []
    for _ in range(REPS):
        st = cs.simulate(N_SIM, lam, eta, 0.0, rng)
        Y.append([2.0 * st["E"] / N_SIM, st["T"] / N_SIM,
                  st["S2"] / N_SIM, st["H3"] / N_SIM, st["H11"] / N_SIM])
    Y = np.array(Y)
    return N_SIM * np.cov(Y.T), Y.mean(0)

def main():
    mm = MomentMap(N=60)
    rng = np.random.default_rng(20260818)
    out = []
    for lam in (1.0, 3.0, 5.0):
        for eta in (5.0, 9.0, 15.0):
            G = G_at(mm, lam, eta)
            Om, ybar = omega_hat(lam, eta, rng)
            I = G.T @ np.linalg.solve(Om, G)
            ieff = 1.0 / np.linalg.inv(I)[2, 2]
            out.append(dict(lam=lam, eta=eta, ieff=float(ieff),
                            lam_max_omega=float(np.linalg.eigvalsh(Om)[-1])))
            print(f"[ieff] lam={lam} eta={eta}: I_eff={ieff:.4f} "
                  f"lam_max(Omega)={out[-1]['lam_max_omega']:.1f}",
                  flush=True)
    json.dump(out, open(os.path.join(BASE, "runs", "ieff_profile.json"),
                        "w"), indent=1)
    print("[ieff] wrote runs/ieff_profile.json", flush=True)

if __name__ == "__main__":
    main()

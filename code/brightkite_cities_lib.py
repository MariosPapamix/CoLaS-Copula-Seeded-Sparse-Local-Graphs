"""Shared per-metro analysis battery used by brightkite_cities.py and
gowalla_cities.py: observed moments, Spearman mixing, degree-preserving
configuration null, direct psi, graph-only fit + bootstrap, compensated
ablation, model-implied channel decomposition, Jacobian diagnostic."""
import sys
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
from scipy.stats import rankdata
from brightkite_pipeline import (region_graph, graph_stats, spatial_score,
                                 CondModel, fit_cond, direct_psi,
                                 sample_marks_cond)

def spearman_r(deg, R, C):
    rk = rankdata(deg, method="average")
    d1 = np.concatenate([rk[R], rk[C]]); d2 = np.concatenate([rk[C], rk[R]])
    return float(np.corrcoef(d1, d2)[0, 1])

def config_null(n, R, C, rng, n_rep=20, swap_mult=15):
    base = np.stack([R, C], 1)
    outs = []
    for _ in range(n_rep):
        e = base.copy()
        eset = set((int(a) << 20) | int(b) for a, b in e)
        m = e.shape[0]; n_sw = swap_mult * m
        idx = rng.integers(0, m, 2 * n_sw)
        done = 0; k = 0
        while done < n_sw and k + 1 < idx.size:
            i, j = idx[k], idx[k + 1]; k += 2
            a, b = e[i]; c, d = e[j]
            if len({a, b, c, d}) < 4: continue
            na, nb = (min(a, d), max(a, d)), (min(c, b), max(c, b))
            k1 = (na[0] << 20) | na[1]; k2 = (nb[0] << 20) | nb[1]
            if k1 in eset or k2 in eset: continue
            eset.discard((min(a, b) << 20) | max(a, b))
            eset.discard((min(c, d) << 20) | max(c, d))
            eset.add(k1); eset.add(k2)
            e[i] = na; e[j] = nb; done += 1
        st = graph_stats(n, e[:, 0], e[:, 1])
        outs.append([st["C"], st["r"]])
    outs = np.array(outs)
    return outs.mean(0), outs.std(0)

def model_channels(model, theta, rng, n_rep=30):
    import scipy.sparse as sp
    lam, Rkm, b = theta
    within = model.pd <= Rkm
    P0, P1 = model.P[within, 0], model.P[within, 1]
    n = model.n
    acc = []
    for _ in range(n_rep):
        V = sample_marks_cond(model.phi, b, rng)
        p = -np.expm1(-lam * V[P0] * V[P1])
        Pm = sp.csr_matrix((np.concatenate([p, p]),
                            (np.concatenate([P0, P1]),
                             np.concatenate([P1, P0]))), shape=(n, n))
        Lam = np.asarray(Pm.sum(1)).ravel()
        keep = rng.random(p.size) < p
        R, C = P0[keep], P1[keep]
        if R.size < 10: continue
        G2 = Pm @ Pm
        Gam = np.asarray(G2[R, C]).ravel()
        L1 = np.concatenate([Lam[R], Lam[C]]); L2 = np.concatenate([Lam[C], Lam[R]])
        D = L1.var() + L1.mean()
        acc.append([np.cov(L1, L2)[0, 1] / D, np.concatenate([Gam, Gam]).mean() / D])
    acc = np.array(acc)
    return acc.mean(0), acc.std(0)

def jac_diag(model, theta, n_rep=40):
    def m(th):
        return model.sim(th[0], th[1], max(th[2], 0.0),
                         np.random.default_rng(4242), n_rep).mean(0)
    th = np.array(theta, float)
    J = np.zeros((3, 3))
    for j in range(3):
        e = np.zeros(3); e[j] = np.array([0.05, 0.05, 0.02])[j]
        J[:, j] = (m(th + e) - m(th - e)) / (2 * e[j])
    sv = np.linalg.svd(J, compute_uv=False)
    return dict(sv=[float(x) for x in sv], smin=float(sv[-1]),
                cond=float(sv[0] / max(sv[-1], 1e-12)))

def run_city(name, users, lats, lons, counts, E, rng, boot=30):
    g = region_graph(name, users, lats, lons, counts, E)
    if g["n"] < 300 or g["R"].size < 200:
        return None
    st = graph_stats(g["n"], g["R"], g["C"])
    deg = np.bincount(g["R"], minlength=g["n"]) + np.bincount(g["C"], minlength=g["n"])
    r_s = spearman_r(deg, g["R"], g["C"])
    (Cn, rn), (Cns, rns) = config_null(g["n"], g["R"], g["C"], rng)
    phi, bw = spatial_score(g["X"])
    b_h, se_b, psi_d, se_pd = direct_psi(g["V"], phi)
    model = CondModel(g["X"], phi)
    target = np.array([st["ell"], st["C"], st["r"]])
    th, q = fit_cond(model, target, np.array([1.0, 1.0, max(b_h, 0.02)]),
                     np.random.default_rng(1))
    th0, _ = fit_cond(model, target, th, np.random.default_rng(2), fix_b=0.0)
    m_fit = model.sim(*th, np.random.default_rng(3), 60).mean(0)
    m_geo = model.sim(*th0, np.random.default_rng(3), 60).mean(0)
    psis = []
    for _ in range(boot):
        mb = model.sim(*th, rng, 1)[0]
        thb, _ = fit_cond(model, mb, th, rng, n_rep=10)
        psis.append(thb[2] ** 2)
    psis = np.array(psis)
    ch_mean, ch_sd = model_channels(model, th, rng)
    jd = jac_diag(model, th)
    print(f"{name}: n={st['n']} ell={st['ell']:.2f} C={st['C']:.3f} "
          f"r={st['r']:+.3f} r_sp={r_s:+.3f} | null C={Cn:.3f} r={rn:+.3f} | "
          f"psi_d={psi_d:.4f}({se_pd:.4f}) psi_g={th[2]**2:.4f} "
          f"[{np.quantile(psis,0.025):.3f},{np.quantile(psis,0.975):.3f}] | "
          f"r_geo={m_geo[2]:+.3f} smin={jd['smin']:.4f}", flush=True)
    return dict(n=st["n"], edges=st["edges"], ell=st["ell"], C=st["C"],
                r=st["r"], r_spearman=r_s,
                config_null=dict(C=float(Cn), C_sd=float(Cns),
                                 r=float(rn), r_sd=float(rns)),
                b_direct=float(b_h), se_b=float(se_b),
                psi_direct=float(psi_d), se_psi=float(se_pd),
                theta=[float(x) for x in th], psi_graph=float(th[2] ** 2),
                psi_boot_q=[float(np.quantile(psis, x)) for x in (0.025, 0.975)],
                model_moments=[float(x) for x in m_fit],
                r_geo=float(m_geo[2]),
                channels=dict(sort=float(ch_mean[0]), overlap=float(ch_mean[1])),
                jac=jd, bw_km=float(bw))

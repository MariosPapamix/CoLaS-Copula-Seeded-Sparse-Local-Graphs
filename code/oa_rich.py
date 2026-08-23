"""Richer-class fit on the OpenAlex Netherlands subsample (two-view
analysis, graph side).  Same protocol as rich_fit_city.py: eight-moment
simulated moment matching, parametric-bootstrap GOF before attribution,
within-class counterfactuals, signed-b profile.

Usage: python3 oa_rich.py free      (b free, signed)
       python3 oa_rich.py fixdir    (b pinned to the direct estimate)
env RICH_SWEEPS / RICH_REPS / RICH_NM as in rich_fit_city.py.
Writes runs/oa_rich_NLsub_<mode>.json
"""
import os, sys, json, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
from scipy.stats import rankdata
import brightkite_pipeline as bp
from richer_model import RichModel, TARGET_KEYS, fit_rich, gof
from oa_analysis import build_nlsub, phi_floor

MODE = sys.argv[1] if len(sys.argv) > 1 else "free"
_BASE = os.environ.get("COLAS_BASE") or os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))

def main():
    g, _ = build_nlsub()
    st = bp.graph_stats(g["n"], g["R"], g["C"])
    phi, bw = phi_floor(g["X"])
    b_dir, se_b, psi_d, se_pd = bp.direct_psi(g["V"], phi)
    d = np.linalg.norm(g["X"][g["R"]] - g["X"][g["C"]], axis=1)
    deg = (np.bincount(g["R"], minlength=g["n"])
           + np.bincount(g["C"], minlength=g["n"]))
    q = np.percentile(d, [25, 50, 75, 90])
    target = [st["ell"], st["C"], st["r"], q[0], q[1], q[2],
              float(deg.std() / deg.mean()), float(np.percentile(deg, 99))]
    rk = rankdata(deg, method="average")
    d1 = np.concatenate([rk[g["R"]], rk[g["C"]]])
    d2 = np.concatenate([rk[g["C"]], rk[g["R"]]])
    obs_held = dict(q90=float(q[3]),
                    hub_share=float(np.sort(deg)[::-1][:max(1, g["n"] // 100)]
                                    .sum() / max(deg.sum(), 1)),
                    max_deg=int(deg.max()),
                    spearman=float(np.corrcoef(d1, d2)[0, 1]))
    print(f"[rich-oa {MODE}] n={g['n']} b_dir={b_dir:+.4f} targets="
          f"{[round(t, 3) for t in target]}", flush=True)

    M = RichModel(g["X"], phi)
    # data-driven start: R from q25, L from q75
    x0 = (2.0, float(np.clip(q[0], 0.3, 3.0)), 0.05,
          float(np.clip(q[2] / 2.0, 5.0, 60.0)), 1.0, 0.05,
          float(b_dir) if MODE == "fixdir" else 0.0)
    if len(sys.argv) > 2:            # alternative start, e.g. mode "free2"
        x0 = tuple(json.loads(sys.argv[2]))
    t0 = time.time()
    th = fit_rich(M, target,
                  sweeps=int(os.environ.get("RICH_SWEEPS", "3")),
                  n_rep=int(os.environ.get("RICH_REPS", "4")), x0=x0,
                  fix_b=(float(b_dir) if MODE == "fixdir" else None))
    print(f"[rich-oa {MODE}] fit done in {time.time()-t0:.0f}s", flush=True)

    G = gof(M, th, target, B=100)
    print(f"[rich-oa {MODE}] GOF max|z| = {G['max_abs_z']:.2f} z="
          f"{[round(z, 1) for z in G['z']]}", flush=True)

    rng = np.random.default_rng(7)
    fulls = M.sim(tuple(th), rng, 20, full=True)
    held = dict(hub_share=[], max_deg=[], spearman=[])
    for f in fulls:
        dg = f["deg"]
        ds = np.sort(dg)[::-1]
        held["hub_share"].append(float(ds[:max(1, M.n // 100)].sum()
                                       / max(dg.sum(), 1)))
        held["max_deg"].append(float(dg.max()))
        A = f["A"].tocoo()
        up = A.row < A.col
        rr, cc = A.row[up], A.col[up]
        rkm = rankdata(dg, method="average")
        e1 = np.concatenate([rkm[rr], rkm[cc]])
        e2 = np.concatenate([rkm[cc], rkm[rr]])
        held["spearman"].append(float(np.corrcoef(e1, e2)[0, 1])
                                if e1.std() > 0 else 0.0)
    held_summary = {k: [float(np.mean(v)), float(np.std(v))]
                    for k, v in held.items() if v}

    def recal_lam(theta):
        theta = list(theta)
        for _ in range(8):
            m = M.sim(tuple(theta), np.random.default_rng(999), 6)
            f = target[0] / max(m.mean(0)[0], 1e-9)
            if abs(f - 1) < 0.01:
                break
            theta[0] = float(np.clip(theta[0] * f, 1e-3, 500))
        return theta

    def r_at(theta, reps=40):
        m = M.sim(tuple(theta), np.random.default_rng(31), reps)
        return float(m.mean(0)[2]), float(m.std(0)[2] / np.sqrt(len(m)))

    r_hat, r_hat_se = r_at(th)
    lam, R, w, L, sigma, tau, b = th
    cf = {}
    for name, theta_cf in [
            ("no_closure",    (lam, R, w, L, sigma, 0.0, b)),
            ("no_hubs",       (lam, R, w, L, 0.05, tau, b)),
            ("no_longrange",  (lam, R, 0.0, L, sigma, tau, b)),
            ("b_zero",        (lam, R, w, L, sigma, tau, 0.0)),
            ("b_direct",      (lam, R, w, L, sigma, tau, float(b_dir))),
            ("b_direct_flip", (lam, R, w, L, sigma, tau, -float(b_dir)))]:
        tcf = recal_lam(theta_cf)
        r_cf, r_cf_se = r_at(tcf)
        cf[name] = dict(r=r_cf, se=r_cf_se, delta=r_hat - r_cf,
                        theta=[float(x) for x in tcf])
        print(f"[rich-oa {MODE}] cf {name}: r={r_cf:+.4f} "
              f"delta={r_hat-r_cf:+.4f}", flush=True)

    prof = []
    for bb in np.linspace(-0.3, 0.3, 9):
        tcf = recal_lam((lam, R, w, L, sigma, tau, float(bb)))
        r_b, se_bp = r_at(tcf, reps=25)
        prof.append([float(bb), r_b, se_bp])
        print(f"[rich-oa {MODE}] profile b={bb:+.3f} r={r_b:+.4f}",
              flush=True)

    out = dict(region="NL-sub", mode=MODE, n=int(g["n"]),
               target=dict(zip(TARGET_KEYS, [float(t) for t in target])),
               obs_held=obs_held,
               theta=dict(zip(["lam", "R", "w", "L", "sigma", "tau", "b"],
                              [float(x) for x in th])),
               b_direct=float(b_dir), se_b_direct=float(se_b),
               gof=G, held=held_summary, r_hat=r_hat, r_hat_se=r_hat_se,
               counterfactuals=cf, b_profile=prof, bw_km=float(bw))
    path = os.path.join(_BASE, "runs", f"oa_rich_NLsub_{MODE}.json")
    json.dump(out, open(path, "w"), indent=1)
    print(f"[rich-oa {MODE}] wrote {path}", flush=True)

if __name__ == "__main__":
    main()

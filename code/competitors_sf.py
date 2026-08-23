"""Competitor suite at SF Bay (Brightkite): DC-SBM, spatial logistic
with node effects, and the closure-augmented richer class, evaluated on
(a) held-out edge prediction (80/20 dyad split, AUC) and
(b) posterior-predictive motif summaries (ell, C, r, edge-length
quantiles, degree CV/p99) against the observed graph.
Writes runs/competitors_sf.json.

Models:
  1. DC-SBM (Karrer--Newman Poisson degree-corrected), K chosen from
     {2,...,8} by held-out AUC; blocks from regularized spectral
     clustering on the training graph, theta/omega by closed-form MLE.
  2. Spatial logistic with node effects:
     logit p_ij = a_i + a_j + beta_0 - beta_1 d_ij, fitted by MAP with
     ridge on a (case-control sampling of non-edges), full-dyad eval.
  3. Richer class of richer_model.py at the fitted theta (from
     runs/rich_brightkite_SFBay.json); simulation-based predictive
     p_ij estimated from repeated realizations for AUC.
  4. Geometry-only baseline (hard disc at fitted (lam,R), psi=0)
     for reference.
"""
import os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import scipy.sparse as sp
import brightkite_pipeline as bp
from scipy.spatial.distance import pdist, squareform

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
rng = np.random.default_rng(20260817)

def auc(scores_pos, scores_neg):
    x = np.concatenate([scores_pos, scores_neg])
    y = np.concatenate([np.ones(scores_pos.size), np.zeros(scores_neg.size)])
    o = np.argsort(x, kind="stable")
    ranks = np.empty_like(o, dtype=float); ranks[o] = np.arange(1, x.size + 1)
    # midranks for ties
    xs = x[o]
    i = 0
    while i < xs.size:
        j = i
        while j + 1 < xs.size and xs[j + 1] == xs[i]:
            j += 1
        if j > i:
            ranks[o[i:j + 1]] = (i + j + 2) / 2.0
        i = j + 1
    rp = ranks[y == 1]
    n1, n0 = scores_pos.size, scores_neg.size
    return float((rp.sum() - n1 * (n1 + 1) / 2) / (n1 * n0))

def graph_stats_edges(n, R_idx, C_idx, X):
    deg = np.bincount(R_idx, minlength=n) + np.bincount(C_idx, minlength=n)
    A = sp.csr_matrix((np.ones(2 * R_idx.size, dtype=np.float32),
                       (np.concatenate([R_idx, C_idx]),
                        np.concatenate([C_idx, R_idx]))), shape=(n, n))
    S2 = 0.5 * float((deg.astype(float) * (deg - 1)).sum())
    T = float((A @ A).multiply(A).sum()) / 6.0
    d1 = np.concatenate([deg[R_idx], deg[C_idx]]).astype(float)
    d2 = np.concatenate([deg[C_idx], deg[R_idx]]).astype(float)
    r = float(np.corrcoef(d1, d2)[0, 1]) if d1.std() > 0 else 0.0
    dd = np.linalg.norm(X[R_idx] - X[C_idx], axis=1)
    q = np.percentile(dd, [25, 50, 75]) if dd.size else [0, 0, 0]
    return dict(ell=2 * R_idx.size / n, C=3 * T / S2 if S2 > 0 else 0,
                r=r, q25=float(q[0]), q50=float(q[1]), q75=float(q[2]),
                deg_cv=float(deg.std() / max(deg.mean(), 1e-9)),
                p99=float(np.percentile(deg, 99)))

def spectral_blocks(A, K, n):
    """Regularized spectral clustering."""
    deg = np.asarray(A.sum(1)).ravel()
    tau = deg.mean()
    Dm = sp.diags(1.0 / np.sqrt(deg + tau))
    Lr = Dm @ A @ Dm
    from scipy.sparse.linalg import eigsh
    k_eig = min(K, n - 2)
    vals, vecs = eigsh(Lr, k=k_eig, which="LA")
    Xs = vecs / np.maximum(np.linalg.norm(vecs, axis=1, keepdims=True), 1e-12)
    from scipy.cluster.vq import kmeans2
    _, lab = kmeans2(Xs, K, minit="++", seed=7)
    return lab

def dcsbm_fit_predict(n, Rtr, Ctr, K):
    A = sp.csr_matrix((np.ones(2 * Rtr.size), (np.concatenate([Rtr, Ctr]),
                       np.concatenate([Ctr, Rtr]))), shape=(n, n))
    lab = spectral_blocks(A, K, n)
    deg = np.asarray(A.sum(1)).ravel()
    m_rs = np.zeros((K, K))
    for a, b in zip(Rtr, Ctr):
        m_rs[lab[a], lab[b]] += 1
        m_rs[lab[b], lab[a]] += 1
    kap = np.array([deg[lab == r].sum() for r in range(K)])
    th = deg / np.maximum(kap[lab], 1e-12)      # theta normalized per block
    def pair_rate(i, j):
        return th[i] * th[j] * m_rs[lab[i], lab[j]]
    return lab, th, m_rs, pair_rate

def main():
    users, lats, lons, counts, E = bp.load_homes_csv(
        os.path.join(BASE, "data/loc-brightkite_edges.txt.gz"),
        os.path.join(BASE, "data/brightkite_homes.csv.gz"))
    bp.REGIONS.update({"SF Bay": (37.15, 38.10, -122.75, -121.60)})
    g = bp.region_graph("SF Bay", users, lats, lons, counts, E)
    n = g["n"]; X = g["X"]
    E_idx = np.stack([g["R"], g["C"]], 1)
    m = E_idx.shape[0]
    D = squareform(pdist(X)).astype(np.float32)

    # 80/20 edge split + matched non-edge test sample
    perm = rng.permutation(m)
    ntest = m // 5
    test_e = E_idx[perm[:ntest]]
    train_e = E_idx[perm[ntest:]]
    eset = set((int(a) << 21) | int(b) for a, b in E_idx)
    neg = []
    while len(neg) < 10 * ntest:
        a = rng.integers(0, n); b = rng.integers(0, n)
        if a == b: continue
        a, b = (a, b) if a < b else (b, a)
        if ((int(a) << 21) | int(b)) not in eset:
            neg.append((a, b))
    neg = np.array(neg)
    Rtr, Ctr = train_e[:, 0], train_e[:, 1]
    res = {}

    # ---- 1. DC-SBM over K ------------------------------------------------
    best = None
    for K in (2, 3, 4, 5, 6, 8):
        lab, th, m_rs, pr = dcsbm_fit_predict(n, Rtr, Ctr, K)
        sp_ = np.array([pr(a, b) for a, b in test_e])
        sn_ = np.array([pr(a, b) for a, b in neg])
        a_ = auc(sp_, sn_)
        if best is None or a_ > best[1]:
            best = (K, a_, lab, th, m_rs)
    K, auc_dcsbm, lab, th, m_rs = best
    # posterior predictive draws
    pp = []
    kap = np.array([np.asarray(th)[lab == r_].sum() for r_ in range(K)])
    for _ in range(25):
        lam_full = np.outer(th, th) * m_rs[np.ix_(lab, lab)]
        U = np.triu_indices(n, 1)
        pe = 1 - np.exp(-lam_full[U])
        keep = rng.random(pe.size) < pe
        pp.append(graph_stats_edges(n, U[0][keep], U[1][keep], X))
    res["dcsbm"] = dict(K=K, auc=auc_dcsbm,
                        pp={k: [float(np.mean([q[k] for q in pp])),
                                float(np.std([q[k] for q in pp]))]
                            for k in pp[0]})
    print("[comp] DC-SBM", res["dcsbm"]["K"], "auc", round(auc_dcsbm, 3),
          flush=True)

    # ---- 2. spatial logistic with node effects (sklearn, case-control) --
    from sklearn.linear_model import LogisticRegression
    ncc = 40 * Rtr.size
    cc = []
    while len(cc) < ncc:
        a = rng.integers(0, n); b = rng.integers(0, n)
        if a == b: continue
        a, b = (a, b) if a < b else (b, a)
        if ((int(a) << 21) | int(b)) not in eset:
            cc.append((a, b))
    cc = np.array(cc)
    pairs = np.concatenate([train_e, cc])
    y = np.concatenate([np.ones(Rtr.size), np.zeros(cc.shape[0])])
    dpair = D[pairs[:, 0], pairs[:, 1]].astype(float)
    rows = np.repeat(np.arange(pairs.shape[0]), 2)
    cols = pairs.ravel()
    inc = sp.csr_matrix((np.ones(rows.size), (rows, cols)),
                        shape=(pairs.shape[0], n))
    Xd = sp.hstack([inc, sp.csr_matrix(dpair[:, None]),
                    sp.csr_matrix(np.log1p(dpair)[:, None])]).tocsr()
    clf = LogisticRegression(C=1.0, solver="liblinear", max_iter=200)
    clf.fit(Xd, y)
    a_eff = clf.coef_[0][:n]
    b_d, b_ld = clf.coef_[0][n], clf.coef_[0][n + 1]
    b0 = clf.intercept_[0]
    def logit_score(a, b):
        d_ = D[a, b]
        return a_eff[a] + a_eff[b] + b0 + b_d * d_ + b_ld * np.log1p(d_)
    sp_ = np.array([logit_score(a, b) for a, b in test_e])
    sn_ = np.array([logit_score(a, b) for a, b in neg])
    auc_log = auc(sp_, sn_)
    U = np.triu_indices(n, 1)
    dU_ = D[U].astype(float)
    etaf = a_eff[U[0]] + a_eff[U[1]] + b0 + b_d * dU_ + b_ld * np.log1p(dU_)
    from scipy.optimize import brentq
    def f(c):
        return (1 / (1 + np.exp(-(etaf + c)))).sum() - Rtr.size / 0.8
    try:
        c0 = brentq(f, -15, 15)
    except ValueError:
        c0 = 0.0
    pf = 1 / (1 + np.exp(-(etaf + c0)))
    pp = []
    for _ in range(25):
        keep = rng.random(pf.size) < pf
        pp.append(graph_stats_edges(n, U[0][keep], U[1][keep], X))
    res["spatial_logistic"] = dict(auc=auc_log, beta_km=float(-b_d),
                                   pp={k: [float(np.mean([q[k] for q in pp])),
                                           float(np.std([q[k] for q in pp]))]
                                       for k in pp[0]})
    print("[comp] logistic auc", round(auc_log, 3), flush=True)

    # ---- 3. richer class (predictive by simulation) ----------------------
    rich_path = os.path.join(BASE, "runs", "rich_brightkite_SFBay.json")
    if os.path.exists(rich_path):
        from richer_model import RichModel
        rich = json.load(open(rich_path))
        thr = tuple(rich["theta"][k] for k in
                    ("lam", "R", "w", "L", "sigma", "tau", "b"))
        phi, bw = bp.spatial_score(X)
        M = RichModel(X, phi)
        cnt = np.zeros((n, n), dtype=np.float32)
        NS = 60
        for _ in range(NS):
            f_ = M.sim(thr, rng, 1, full=True)
            if not f_:
                continue
            A_ = f_[0]["A"]
            cnt += A_.toarray()
        pfreq = cnt / NS
        sp_ = np.array([pfreq[a, b] for a, b in test_e])
        sn_ = np.array([pfreq[a, b] for a, b in neg])
        auc_rich = auc(sp_, sn_)
        # marks-conditional predictive score: observed V + fitted kernel
        Vobs = g["V"]
        kapf = (D <= thr[1]).astype(float) + thr[2] * np.exp(-D / max(thr[3], 1e-3))
        def rich_cond(a, b):
            return Vobs[a] * Vobs[b] * kapf[a, b]
        spc = np.array([rich_cond(a, b) for a, b in test_e])
        snc = np.array([rich_cond(a, b) for a, b in neg])
        auc_rich_cond = auc(spc, snc)
        pp = []
        for _ in range(25):
            f_ = M.sim(thr, rng, 1, full=True)
            if f_:
                A_ = f_[0]["A"].tocoo(); up = A_.row < A_.col
                pp.append(graph_stats_edges(n, A_.row[up], A_.col[up], X))
        res["richer_class"] = dict(auc=auc_rich,
                                   auc_marks_conditional=auc_rich_cond,
                                   pp={k: [float(np.mean([q[k] for q in pp])),
                                           float(np.std([q[k] for q in pp]))]
                                       for k in pp[0]})
        print("[comp] richer auc", round(auc_rich, 3), flush=True)

    # ---- 4. geometry-only baseline ---------------------------------
    r8 = json.load(open(os.path.join(BASE, "runs",
                                     "cities_r8_brightkite.json")))
    th8 = r8["Brightkite/SF Bay"]["theta"]
    lam8, R8 = th8[0], th8[1]
    dU = D[U]
    within = dU <= R8
    Vm = rng.random(n)  # marks resampled per draw below
    pp = []
    prob_cache = None
    for _ in range(25):
        V = rng.random(n)
        pe = np.where(within, 1 - np.exp(-lam8 * V[U[0]] * V[U[1]]), 0.0)
        keep = rng.random(pe.size) < pe
        pp.append(graph_stats_edges(n, U[0][keep], U[1][keep], X))
    pe_mean = np.where(within, 1 - np.exp(-lam8 * 0.25), 0.0)  # E[VV']=1/4
    sp_ = np.array([pe_mean[np.searchsorted(np.ravel_multi_index((U[0], U[1]), (n, n)),
                    np.ravel_multi_index((a, b), (n, n)))] if False else
                    (1 - np.exp(-lam8 * 0.25)) * (D[a, b] <= R8)
                    for a, b in test_e])
    sn_ = np.array([(1 - np.exp(-lam8 * 0.25)) * (D[a, b] <= R8)
                    for a, b in neg])
    res["colas_geo"] = dict(auc=auc(sp_, sn_),
                            pp={k: [float(np.mean([q[k] for q in pp])),
                                    float(np.std([q[k] for q in pp]))]
                                for k in pp[0]})
    print("[comp] colas_geo auc", round(res["colas_geo"]["auc"], 3), flush=True)

    obs = graph_stats_edges(n, g["R"], g["C"], X)
    res["observed"] = {k: float(v) for k, v in obs.items()}
    res["split"] = dict(train_edges=int(Rtr.size), test_edges=int(ntest),
                        neg_samples=int(neg.shape[0]))
    json.dump(res, open(os.path.join(BASE, "runs", "competitors_sf.json"),
                        "w"), indent=1)
    print("[comp] wrote runs/competitors_sf.json", flush=True)

if __name__ == "__main__":
    main()

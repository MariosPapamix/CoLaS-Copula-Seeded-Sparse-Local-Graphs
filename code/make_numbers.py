"""Generate revision/numbers.tex from simulation + empirical results."""
import sys, json, os
import os as _os
_BASE = _os.environ.get("COLAS_BASE") or _os.path.dirname(
    _os.path.dirname(_os.path.abspath(__file__)))
def _p(*parts):
    q = _os.path.join(_BASE, *parts)
    _os.makedirs(_os.path.dirname(q), exist_ok=True) if parts[:1] in (("runs",), ("figures",), ("revision",)) else None
    return q
import numpy as np

RUNS = _p("runs") + _os.sep
OUT = _p("revision", "numbers.tex")
R = np.load(RUNS + "sim_results.npz", allow_pickle=True)

def pct(x, d=1):
    return f"{100*x:.{d}f}"

L = []
def cmd(name, val):
    """Define or overwrite a macro (later sections may refresh earlier ones)."""
    tag = f"\\newcommand{{\\{name}}}{{"
    line = f"\\newcommand{{\\{name}}}{{{val}}}"
    for _i, _existing in enumerate(L):
        if _existing.startswith(tag):
            L[_i] = line
            return
    L.append(line)

# coverage table
for n, tag in ((1000, "Thousand"), (10000, "Tenk"), (100000, "Hundk")):
    hit = R[f"cov_hit_{n}"]; th = R[f"cov_theta_{n}"]
    bias = th.mean(0) - np.array([3.0, 9.0, 0.0835])
    for j, p in enumerate(("Lam", "Eta", "Psi")):
        cmd(f"cov{p}{tag}", pct(hit[:, j].mean()))
        b = bias[j]
        cmd(f"bias{p}{tag}", f"{b:+.4f}" if abs(b) >= 5e-5 else f"{b:+.5f}")
cmd("repSmall", "400"); cmd("repLarge", "100")
cmd("covPsiSmallPct", pct(R["cov_hit_1000"][:, 2].mean()) + "\\%")

# moment bias (r-coordinate)
cmd("biasRthousand", f"{1000*R['moment_bias_1000'][2]:.0f}")
cmd("biasRtenk", f"{10000*R['moment_bias_10000'][2]:.0f}")

# feasible arm
cmd("feasCov", " / ".join(pct(x) for x in R["feas_hit"].mean(0)))

# boundary
D = R["null_D_10000"]
cmd("nullMassTenkPct", pct(np.mean(D <= 1e-12)) + "\\%")
cmd("sizeTenkPct", pct(np.mean(D > 2.705543)) + "\\%")
cmd("Ieff", f"{float(R['I_eff_null']):.3f}")

# power
psis = R["power_psis"]; pw = R["power"]
for ps, name in ((0.02, "powTwoPct"), (0.04, "powFourPct"), (0.0835, "powFullPct")):
    i = int(np.argmin(np.abs(psis - ps)))
    cmd(name, pct(pw[i]) + "\\%")

# misspecification
tri = R["miss_tri"].mean(0); bump = R["miss_bump"].mean(0); twoh = R["miss_twoharm"].mean(0)
cmd("missKernelSummary",
    f"$-{100*(0.0835-tri[2])/0.0835:.0f}\\%$ (triangular) and "
    f"$-{100*(0.0835-bump[2])/0.0835:.0f}\\%$ (truncated Gaussian)")
cmd("missTwoharmPsi", f"$\\widehat\\psi={twoh[2]:.4f}$")

# ---- demo benchmark (always available after run) ----------------------
dm = RUNS + "demo_results.json"
if os.path.exists(dm):
    with open(dm) as fh:
        Dm = json.load(fh)
    fitn = Dm["fit_region"]; heldn = Dm["held_region"]
    st = Dm[fitn]["stats"]
    cmd("dmN", f"{st['n']:,}")
    cmd("dmR", f"{st['r']:.3f}")
    th = Dm["fit"]["theta"]
    cmd("dmLamHat", f"{th[0]:.2f}")
    cmd("dmRHat", f"{th[1]:.2f}")
    cmd("dmPsiGraph", f"{Dm['fit']['psi_graph']:.3f}")
    qlo, qhi = Dm["fit"]["psi_boot_q"]
    cmd("dmPsiCI", f"[{qlo:.3f},{qhi:.3f}]")
    cmd("dmPsiDirect", f"{Dm[fitn]['psi_direct']:.4f}")
    cmd("dmPsiDirectSE", f"{Dm[fitn]['se_psi']:.4f}")
    cmd("dmRgeo", f"{Dm['ablation']['model_moments'][2]:.3f}")
    cmd("dmShare", f"{Dm['ablation']['geometric_share_pct']:.0f}\\%")
    tq = Dm.get("transfer", {})
    obs = np.array(tq["obs"]); mod = np.array(tq["model"])
    rel = float(np.max(np.abs(mod - obs) / np.maximum(np.abs(obs), 1e-9)))
    cmd("dmTransferPct", f"{100*rel:.0f}\\%")
else:
    for nm in ("dmN dmR dmLamHat dmRHat dmPsiGraph dmPsiCI "
               "dmPsiDirect dmPsiDirectSE dmRgeo dmShare dmTransferPct").split():
        cmd(nm, "---")

# ---- Brightkite real fit ----------------------------------------------
bk = RUNS + "brightkite_results.json"
if os.path.exists(bk):
    with open(bk) as fh:
        B = json.load(fh)
    fit_name = B["fit_region"]; held = B["held_region"]
    st = B[fit_name]["stats"]; sh = B[held]["stats"]
    cmd("bkC", f"{st['C']:.3f}")
    cmd("bkR", f"{st['r']:.3f}")
    cmd("bkBDirect", f"{B[fit_name]['b_direct']:.3f}")
    cmd("bkBDirectSE", f"{B[fit_name]['se_b']:.3f}")
    qlo, qhi = B["fit"]["psi_boot_q"]
    cmd("bkRgeo", f"{B['ablation']['model_moments'][2]:.3f}")
    gap = B["ablation"]["model_moments"][2] - st["r"]
    cmd("bkGap", f"{gap:.2f}")
    g = B.get("geo", {}).get(fit_name, {})
    cmd("bkMedKm", f"{g.get('med_km', float('nan')):.1f}")
    cmd("bkPctOneKm", f"{g.get('pct_within_1km', float('nan')):.1f}\\%")
    cmd("bkMaxDeg", f"{g.get('max_deg', 0):,}")
    cmd("bkHubShare", f"{g.get('top1pct_edge_share', float('nan')):.0f}\\%")
else:
    for nm in ("bkC bkR bkBDirect bkBDirectSE bkRgeo bkGap bkMedKm bkPctOneKm "
               "bkMaxDeg bkHubShare").split():
        cmd(nm, "---")

with open(OUT, "w") as fh:
    fh.write("% auto-generated from the replication package -- do not edit\n")
    fh.write("\n".join(L) + "\n")
print("wrote", OUT, f"({len(L)} macros)")
print("\n".join(L[:12]))

# ---- grid, five-moment GMM, boundary nuisance --------------------------
import sys as _sys
import os as _os0
_sys.path.insert(0, _os0.path.dirname(_os0.path.abspath(__file__)))

from colas_moments import MomentMap as _MM, channels as _channels
_mm = _MM(N=60)
R3 = np.load(RUNS + "sim_r3.npz", allow_pickle=True)
with open(RUNS + "sim_r3_extra.json") as _fh:
    X3 = json.load(_fh)

ec = R3["ellipse_cov"]; jn = R3["joint_naive"]
for tag, i in (("Thousand", 0), ("Tenk", 1), ("Hundk", 2)):
    cmd(f"ellipseCov{tag}", pct(ec[i]) + "\\%")
    cmd(f"jointNaive{tag}", pct(jn[i]) + "\\%")

rs0, ro0 = _channels(_mm, 3.0, 9.0, 0.0)
rsC, roC = _channels(_mm, 3.0, 9.0, 0.0835)
cmd("chanSortNull", f"{rs0:.3f}")
cmd("chanSortCenter", f"{rsC:.3f}")
cmd("chanOverlapCenter", f"{roC:.3f}")

# grid table
rows = []
for g in X3["grid"]:
    th = g["theta"]
    rows.append("$({:.1f},{:.0f},{:.3f})$ & {:+.2f} & {:.4f} & {:.1f} & {:.1f} & {:.3f} & {:.1f} & {:.0f}\\\\".format(
        th[0], th[1], th[2], g["share0"], g["rmse"][2],
        100*g["cov_marg"][2], 100*g["cov_ell"], g["len_psi"],
        100*g["share_cov"], 100*g["ok"]))
tab = ("\\begin{tabular}{@{}lccccccc@{}}\\toprule\n"
       "$(\\lambda,\\eta,\\psi)$ & $s_0$ & RMSE$_\\psi$ & cov$_\\psi$ & ellipse & len$_\\psi$ & cov$_{s}$ & solve\\% \\\\ \\midrule\n"
       + "\n".join(rows) + "\n\\bottomrule\\end{tabular}")
cmd("gridTable", tab.replace("\n", " "))

th3 = R3["fm_th3"]; th5 = R3["fm_th5"]
t0v = np.array([3.0, 9.0, 0.0835])
cmd("rmseThreePsi", f"{np.sqrt(((th3-t0v)**2).mean(0))[2]:.4f}")
cmd("rmseFivePsi", f"{np.sqrt(((th5-t0v)**2).mean(0))[2]:.4f}")
cmd("covFivePct", "/".join(pct(x) for x in R3["fm_cov5"].mean(0)) + "\\%")
from scipy import stats as _st
cmd("jTestPct", pct(np.mean(R3["fm_J"] > _st.chi2.ppf(0.95, 2))) + "\\%")

shT = R3["pair_shares_true"]
cmd("pairShareA", f"{shT[0]:.2f}")
cmd("pairShareB", f"{shT[1]:.2f}")
cmd("pairShareHatA", f"{R3['pairA_share'].mean():.2f}")
cmd("pairShareHatB", f"{R3['pairB_share'].mean():.2f}")

bn = X3["boundary_nuisance"]
ks = sorted(bn.keys())
cmd("sizeNuisA", pct(bn[ks[0]]["size"]) + "\\%")
cmd("sizeNuisB", pct(bn[ks[1]]["size"]) + "\\%")
cmd("sizeNuisC", pct(bn[ks[2]]["size"]) + "\\%")
cmd("massNuisC", f"{bn[ks[2]]['mass0']:.2f}")

# ---- cities ------------------------------------------------------------
ct = RUNS + "brightkite_cities.json"
if os.path.exists(ct):
    with open(ct) as fh:
        CT = json.load(fh)
    crows = []
    for name in ("SF Bay", "LA", "NYC", "Denver", "London"):
        if name not in CT: continue
        c = CT[name]
        crows.append("{} & {:,} & {:.1f} & {:.3f} & ${:+.3f}$ & ${:+.3f}$ & ${:+.3f}$ & {:.3f} & {:.4f} & [0,{:.3f}] & ${:+.2f}$ & {:.3f}\\\\".format(
            name, c["n"], c["ell"], c["C"], c["r"], c["r_spearman"],
            c["config_null"]["r"], c["config_null"]["C"],
            c["psi_direct"], c["psi_boot_q"][1], c["r_geo"], c["jac"]["smin"]))
    ctab = ("\\begin{tabular}{@{}lccccccccccc@{}}\\toprule\n"
            "Region & $n$ & $\\widehat\\ell$ & $\\widehat C$ & $\\widehat r$ & $\\widehat r_{\\mathrm{rank}}$ & $r_{\\mathrm{cfg}}$ & $C_{\\mathrm{cfg}}$ & $\\widehat\\psi_{\\mathrm{dir}}$ & $\\widehat\\psi_{\\mathrm{gr}}$ CI & $r_{\\mathrm{geo}}$ & $\\sigma_{\\min}$\\\\ \\midrule\n"
            + "\n".join(crows) + "\n\\bottomrule\\end{tabular}")
    cmd("cityTable", ctab.replace("\n", " "))
    sf = CT["SF Bay"]
    cmd("bkSpearman", f"{sf['r_spearman']:+.3f}")
    cmd("bkCfgC", f"{sf['config_null']['C']:.3f}")
    cmd("bkCfgR", f"{sf['config_null']['r']:+.3f}")
    cmd("bkChSort", f"{sf['channels']['sort']:.3f}")
    cmd("bkChOv", f"{sf['channels']['overlap']:.3f}")
    cmd("bkJacSmin", f"{sf['jac']['smin']:.3f}")
    cmd("bkJacCond", f"{sf['jac']['cond']:.0f}")
    _cities = [CT[k] for k in ("SF Bay","LA","NYC","Denver","London") if k in CT]
    cmd("bkRLo", f"{min(c['r'] for c in _cities):+.2f}")
    cmd("bkRHi", f"{max(c['r'] for c in _cities):+.2f}")
    cmd("bkRgeoLo", f"{min(c['r_geo'] for c in _cities):+.2f}")
    cmd("bkRgeoHi", f"{max(c['r_geo'] for c in _cities):+.2f}")
    cmd("bkPsiDirLo", f"{min(c['psi_direct'] for c in _cities):.4f}")
    cmd("bkPsiDirHi", f"{max(c['psi_direct'] for c in _cities):.4f}")
    if "benchmark_n12000" in CT:
        b12 = CT["benchmark_n12000"]
        cmd("dmTwelveCI", f"[{b12['psi_boot_q'][0]:.3f},{b12['psi_boot_q'][1]:.3f}]")
else:
    for nm in ("cityTable bkSpearman bkCfgC bkCfgR bkChSort bkChOv "
               "bkJacSmin bkJacCond dmTwelveCI").split():
        cmd(nm, "---")

# rewrite output with all macros
with open(OUT, "w") as fh:
    fh.write("% auto-generated from the replication package -- do not edit\n")
    fh.write("\n".join(L) + "\n")
print("rewrote", OUT, f"({len(L)} macros)")

# ---- Hausman test, robustness, distance null ---------------------------
R4 = np.load(RUNS + "sim_r4.npz", allow_pickle=True)
cmd("hausSize", pct(float(R4["size_rej"])) + "\\%")
cmd("hausLRtwo", pct(float(R4["lr2_rej"])) + "\\%")
cmd("hausLRfive", pct(float(R4["lr5_rej"])) + "\\%")
cmd("hausHeavy", pct(float(R4["heavy_rej"])) + "\\%")
cmd("lrTwoPsi", f"{R4['lr2_psig'].mean():.4f}")
cmd("lrFivePsi", f"{R4['lr5_psig'].mean():.4f}")
cmd("heavyPsi", f"{R4['heavy_psig'].mean():.4f}")
with open(RUNS + "brightkite_cities.json") as fh:
    _CT2 = json.load(fh)
if "dist_null" in _CT2.get("SF Bay", {}):
    cmd("bkDistC", f"{_CT2['SF Bay']['dist_null']['C']:.3f}")
    cmd("bkDistR", f"{_CT2['SF Bay']['dist_null']['r']:+.3f}")
else:
    cmd("bkDistC", "---"); cmd("bkDistR", "---")

with open(OUT, "w") as fh:
    fh.write("% auto-generated from the replication package -- do not edit\n")
    fh.write("\n".join(L) + "\n")
print("rewrote", OUT, f"({len(L)} macros)")

# ---- Gowalla, competitor, ANND, boundary test, d2 certificate ----------
sfx = RUNS + "sf_extras.json"
if os.path.exists(sfx):
    with open(sfx) as fh:
        SX = json.load(fh)
    ld = SX["latent_distance"]
    cmd("ldBeta", f"{ld['beta']:.3f}")
    cmd("ldC", f"{ld['moments_mean'][1]:.4f}")
    cmd("ldR", f"{ld['moments_mean'][2]:+.3f}")
    cmd("btSFp", f"{SX['boundary_test_SF']['p']:.2f}")
    cmd("btSize", pct(SX["boundary_test_size"]["size"]) + "\\%")
else:
    for nm in "ldBeta ldC ldR btSFp btSize".split():
        cmd(nm, "---")

gw = RUNS + "gowalla_cities.json"
if os.path.exists(gw):
    with open(gw) as fh:
        GW = json.load(fh)
    grows = []
    order = [k for k in ("Austin","Stockholm","SF Bay","Dallas","NYC","London") if k in GW]
    for name in order:
        c = GW[name]
        grows.append("{} & {:,} & {:.1f} & {:.3f} & ${:+.3f}$ & ${:+.3f}$ & ${:+.3f}$ & {:.3f} & {:.4f} & [0,{:.3f}] & ${:+.2f}$ & {:.3f}\\\\".format(
            name, c["n"], c["ell"], c["C"], c["r"], c["r_spearman"],
            c["config_null"]["r"], c["config_null"]["C"],
            c["psi_direct"], c["psi_boot_q"][1], c["r_geo"], c["jac"]["smin"]))
    gtab = ("\\begin{tabular}{@{}lccccccccccc@{}}\\toprule\n"
            "Region & $n$ & $\\widehat\\ell$ & $\\widehat C$ & $\\widehat r$ & $\\widehat r_{\\mathrm{rank}}$ & $r_{\\mathrm{cfg}}$ & $C_{\\mathrm{cfg}}$ & $\\widehat\\psi_{\\mathrm{dir}}$ & $\\widehat\\psi_{\\mathrm{gr}}$ CI & $r_{\\mathrm{geo}}$ & $\\sigma_{\\min}$\\\\ \\midrule\n"
            + "\n".join(grows) + "\n\\bottomrule\\end{tabular}")
    cmd("gwTable", gtab.replace("\n", " "))
    au = GW["Austin"]
    cmd("gwAusBDir", f"{au['b_direct']:.3f}")
    cmd("gwAusBDirSE", f"{au['se_b']:.3f}")
    cmd("gwAusPsiD", f"{au['psi_direct']:.4f}")
    cmd("gwAusPsiDSE", f"{au['se_psi']:.4f}")
    cmd("gwAusR", f"{au['r']:+.3f}")
    cmd("gwAusSpear", f"$+{au['r_spearman']:.3f}$")
    gcities = [GW[k] for k in order]
    cmd("gwRLo", f"${min(c['r'] for c in gcities):+.2f}$")
    cmd("gwRHi", f"${max(c['r'] for c in gcities):+.2f}$")
    cmd("gwRgeoLo", f"$+{min(c['r_geo'] for c in gcities):.2f}$")
    cmd("gwRgeoHi", f"$+{max(c['r_geo'] for c in gcities):.2f}$")
else:
    for nm in ("gwTable gwAusBDir gwAusBDirSE gwAusPsiD gwAusPsiDSE gwAusR "
               "gwAusSpear gwRLo gwRHi gwRgeoLo gwRgeoHi").split():
        cmd(nm, "---")

with open(OUT, "w") as fh:
    fh.write("% auto-generated from the replication package -- do not edit\n")
    fh.write("\n".join(L) + "\n")
print("rewrote", OUT, f"({len(L)} macros)")

# ---- GOF, rails, signed counterfactuals, city tables, pair deltas ------
r8 = {}
for _f in ("cities_r8_brightkite.json", "cities_r8_gowalla.json",
           "cities_r8.json"):
    _pth = RUNS + _f
    if os.path.exists(_pth):
        with open(_pth) as fh:
            r8.update(json.load(fh))

def _round_up_1sig(x):
    """Smallest 1-significant-figure decimal >= x (for honest bounds)."""
    import math
    if x <= 0: return 0.0
    e = math.floor(math.log10(x))
    m = math.ceil(x / 10**e)
    if m == 10: m, e = 1, e + 1
    return m * 10**e

BK_ORDER = ["SF Bay", "LA", "NYC", "Denver", "London"]
GW_ORDER = ["Austin", "Stockholm", "SF Bay", "Dallas", "NYC", "London"]

if r8 and all(f"Brightkite/{c}" in r8 for c in BK_ORDER) \
      and all(f"Gowalla/{c}" in r8 for c in GW_ORDER):
    bk8 = {c: r8[f"Brightkite/{c}"] for c in BK_ORDER}
    gw8 = {c: r8[f"Gowalla/{c}"] for c in GW_ORDER}
    all8 = list(bk8.values()) + list(gw8.values())
    def _calib_ok(c):
        rl = c["rails"]
        return not (rl["b_high"] or rl["R_low"] or rl["R_high"])
    conv8 = [c for c in all8 if _calib_ok(c)]
    bd = {c["boot_draws"] for c in all8}
    cmd("bootDraws", f"{min(bd)}" if len(bd) == 1 else f"{min(bd)}--{max(bd)}")

    def _ctab(cities, order):
        rows = []
        for name in order:
            c = cities[name]
            rows.append("{} & {:,} & {:.1f} & {:.3f} & ${:+.3f}$ & ${:+.3f}$ & ${:+.3f}$ & {:.3f} & {:.4f} & [0,{:.3f}] & ${:+.2f}$ & {:.3f}\\\\".format(
                name, c["n"], c["ell"], c["C"], c["r"], c["r_spearman"],
                c["config_null"]["r"], c["config_null"]["C"],
                c["psi_direct"], c["psi_boot_q"][1], c["r_geo"],
                c["jac"]["smin"]))
        return ("\\begin{tabular}{@{}lccccccccccc@{}}\\toprule\n"
                "Region & $n$ & $\\widehat\\ell$ & $\\widehat C$ & $\\widehat r$ & $\\widehat r_{\\mathrm{rank}}$ & $r_{\\mathrm{cfg}}$ & $C_{\\mathrm{cfg}}$ & $\\widehat\\psi_{\\mathrm{dir}}$ & $\\widehat\\psi_{\\mathrm{gr}}$ CI & $r_{\\mathrm{geo}}$ & $\\sigma_{\\min}$\\\\ \\midrule\n"
                + "\n".join(rows) + "\n\\bottomrule\\end{tabular}")
    cmd("cityTable", _ctab(bk8, BK_ORDER).replace("\n", " "))
    cmd("gwTable", _ctab(gw8, GW_ORDER).replace("\n", " "))

    # GOF table: platform, city, r_mod, r_obs, resid, z, rails, %>R, dsort
    grows = []
    for plat, cities, order in (("BK", bk8, BK_ORDER), ("GW", gw8, GW_ORDER)):
        for name in order:
            c = cities[name]
            rail_bits = []
            if c["rails"]["b_low"]: rail_bits.append("$b{=}0$")
            if c["rails"]["b_high"]: rail_bits.append("$b{=}0.55$")
            if c["rails"]["R_low"]: rail_bits.append("$R{=}0.1$")
            if c["rails"]["R_high"]: rail_bits.append("$R{=}R_{\\max}$")
            rail = ", ".join(rail_bits) if rail_bits else "interior"
            _ds = "${:+.4f}$".format(c["dsort_signed"]) if _calib_ok(c) else "--"
            grows.append("{} {} & ${:+.3f}$ & ${:+.3f}$ & ${:+.3f}$ & ${:+.0f}$ & {} & {:.0f} & {}\\\\".format(
                plat, name, c["model_moments"][2], c["r"], c["r_resid"],
                c["gof_z"], rail, 100*c["frac_outside"], _ds))
    gtab = ("\\begin{tabular}{@{}lccccccc@{}}\\toprule\n"
            "City & $r_{\\mathrm{mod}}$ & $\\widehat r$ & resid & $z$ & rail & $\\%{>}\\widehat R$ & $\\Delta_{\\mathrm{sort}}(\\widehat b_{\\mathrm{dir}})$\\\\ \\midrule\n"
            + "\n".join(grows) + "\n\\bottomrule\\end{tabular}")
    cmd("gofTable", gtab.replace("\n", " "))

    gaps = [c["r_geo"] - c["r"] for c in all8]
    cmd("gapLoAll", f"+{min(gaps):.2f}")
    cmd("gapHiAll", f"+{max(gaps):.2f}")
    cmd("rgeoLoAll", f"+{min(c['r_geo'] for c in all8):.2f}")
    cmd("rgeoHiAll", f"+{max(c['r_geo'] for c in all8):.2f}")
    cmd("robsLoAll", f"{min(c['r'] for c in all8):+.2f}")
    cmd("robsHiAll", f"{max(c['r'] for c in all8):+.2f}")
    fr = [100*c["frac_outside"] for c in all8]
    cmd("fracOutLoAll", f"{min(fr):.0f}")
    cmd("fracOutHiAll", f"{max(fr):.0f}")
    dmax = max(abs(c["dsort_signed"]) for c in conv8)
    cmd("dsortBoundAll", f"{_round_up_1sig(dmax):g}")
    bkfr = [100*c["frac_outside"] for c in bk8.values()]
    gwfr = [100*c["frac_outside"] for c in gw8.values()]
    cmd("bkFracLo", f"{min(bkfr):.0f}"); cmd("bkFracHi", f"{max(bkfr):.0f}")
    cmd("gwFracLo", f"{min(gwfr):.0f}"); cmd("gwFracHi", f"{max(gwfr):.0f}")
    cmd("bkRhatSF", f"{bk8['SF Bay']['theta'][1]:.2f}")
    cmd("bkDsortMax", f"{_round_up_1sig(max(abs(c['dsort_signed']) for c in bk8.values() if _calib_ok(c))):g}")
    cmd("gwAusDsort", f"{gw8['Austin']['dsort_signed']:+.4f}")
    cmd("gwAusDsortFlip", f"{gw8['Austin']['dsort_flip']:+.4f}")

    # refresh SF-derived macros from the r8 run
    sf8 = bk8["SF Bay"]
    cmd("bkRgeo", f"{sf8['r_geo']:.3f}")
    cmd("bkGap", f"{sf8['r_geo'] - sf8['r']:.2f}")
    _nyc8 = bk8["NYC"]
    cmd("bkNycSmin", f"{_nyc8['jac']['smin']:.3f}")
    cmd("bkNycCond", f"{_nyc8['jac']['cond']:,.0f}")
    _sm = sorted(c["jac"]["smin"] for k, c in bk8.items() if k != "NYC")
    cmd("bkSminLo", f"{_sm[0]:.3f}")
    cmd("bkSminHi", f"{_sm[-1]:.3f}")
    cmd("bkSpearman", f"{sf8['r_spearman']:+.3f}")
    cmd("bkCfgC", f"{sf8['config_null']['C']:.3f}")
    cmd("bkCfgR", f"{sf8['config_null']['r']:+.3f}")
    cmd("bkChSort", f"{sf8['channels']['int_']:.3f}")
    cmd("bkChOv", f"{sf8['channels']['overlap']:.3f}")
    cmd("bkJacSmin", f"{sf8['jac']['smin']:.3f}")
    cmd("bkJacCond", f"{sf8['jac']['cond']:.0f}")
    _bk_all = list(bk8.values())
    cmd("bkRLo", f"{min(c['r'] for c in _bk_all):+.2f}")
    cmd("bkRHi", f"{max(c['r'] for c in _bk_all):+.2f}")
    cmd("bkRgeoLo", f"{min(c['r_geo'] for c in _bk_all):+.2f}")
    cmd("bkRgeoHi", f"{max(c['r_geo'] for c in _bk_all):+.2f}")
    cmd("bkPsiDirLo", f"{min(c['psi_direct'] for c in _bk_all):.4f}")
    cmd("bkPsiDirHi", f"{max(c['psi_direct'] for c in _bk_all):.4f}")
    st8 = gw8["Stockholm"]
    cmd("gwStoR", f"{st8['r']:+.3f}")
    cmd("gwStoCfgR", f"{st8['config_null']['r']:+.3f}")
    cmd("gwStoRgeo", f"{st8['r_geo']:+.3f}")
    au8 = gw8["Austin"]
    cmd("gwAusBDir", f"{au8['b_direct']:.3f}")
    cmd("gwAusBDirSE", f"{au8['se_b']:.3f}")
    cmd("gwAusPsiD", f"{au8['psi_direct']:.4f}")
    cmd("gwAusPsiDSE", f"{au8['se_psi']:.4f}")
    cmd("gwAusR", f"{au8['r']:+.3f}")
    cmd("gwAusSpear", f"$+{au8['r_spearman']:.3f}$")
    _gw_all = list(gw8.values())
    cmd("gwRLo", f"${min(c['r'] for c in _gw_all):+.2f}$")
    cmd("gwRHi", f"${max(c['r'] for c in _gw_all):+.2f}$")
    cmd("gwRgeoLo", f"$+{min(c['r_geo'] for c in _gw_all):.2f}$")
    cmd("gwRgeoHi", f"$+{max(c['r_geo'] for c in _gw_all):.2f}$")
else:
    for nm in ("bootDraws gofTable gapLoAll gapHiAll rgeoLoAll rgeoHiAll "
               "gwStoR gwStoCfgR gwStoRgeo "
               "bkNycSmin bkNycCond bkSminLo bkSminHi "
               "robsLoAll robsHiAll fracOutLoAll fracOutHiAll dsortBoundAll "
               "bkFracLo bkFracHi gwFracLo gwFracHi bkRhatSF bkDsortMax "
               "gwAusDsort gwAusDsortFlip").split():
        cmd(nm, "---")

# pair-experiment sorting contributions (compensated increments)
from scipy.optimize import brentq as _brentq
def _delta_sort(mm_, th):
    lam, eta, psi = th
    t = mm_.targets(lam, eta, psi)
    ell0, C0, r1 = float(t["ell"]), float(t["C"]), float(t["r"])
    def _g(lam2):
        return float(mm_.targets(lam2, 1.0, 0.0)["C"]) - C0
    lam2 = _brentq(_g, 0.3, 8.0, xtol=1e-12)
    E2 = float(mm_.targets(lam2, 1.0, 0.0)["motifs"]["E"])
    eta2 = ell0 / (2.0 * E2)
    r0 = float(mm_.targets(lam2, eta2, 0.0)["r"])
    return r1 - r0
_dA = _delta_sort(_mm, (1.6, 11.0, 0.1))
_dB = _delta_sort(_mm, (4.19962088, 6.28932698, 0.01))
cmd("pairDeltaA", f"{_dA:+.3f}")
cmd("pairDeltaB", f"{_dB:+.3f}")

with open(OUT, "w") as fh:
    fh.write("% auto-generated from the replication package -- do not edit\n")
    fh.write("\n".join(L) + "\n")
print("rewrote", OUT, f"({len(L)} macros)")

# ---- richer class, competitors, information profile --------------------
def _load(pth):
    p_ = RUNS + pth
    return json.load(open(p_)) if os.path.exists(p_) else None

rich_sf = _load("rich_brightkite_SFBay.json")
rich_au = _load("rich_gowalla_Austin.json")
comp = _load("competitors_sf.json")
prof = _load("ieff_profile.json")

if prof:
    vals = [e["ieff"] for e in prof]
    cmd("ieffLo", f"{min(vals):.2f}")
    cmd("ieffHi", f"{max(vals):.2f}")
else:
    cmd("ieffLo", "---"); cmd("ieffHi", "---")

if r8:
    cmd("gofZWorstOld", f"{max(abs(c['gof_z']) for c in all8):.0f}")

_RKEYS = ["ell", "C", "r", "q25", "q50", "q75", "deg_cv", "p99"]
def _rich_rows(name, d):
    tg = [d["target"][k] for k in _RKEYS]
    mn = d["gof"]["mean"]; sd = d["gof"]["sd"]; zz = d["gof"]["z"]
    r1 = name + " & obs & " + " & ".join(
        f"${t:+.3f}$" if k == "r" else f"{t:.2f}"
        for k, t in zip(_RKEYS, tg)) + "\\\\"
    r2 = " & model & " + " & ".join(
        (f"${m:+.3f}$" if k == "r" else f"{m:.2f}") + f" ({s:.2f})"
        for k, m, s in zip(_RKEYS, mn, sd)) + "\\\\"
    r3 = " & $z$ & " + " & ".join(f"${z:+.1f}$" for z in zz) + "\\\\"
    return [r1, r2, r3]

if rich_sf and rich_au:
    rows = _rich_rows("SF Bay", rich_sf) + ["\\midrule"] +            _rich_rows("Austin", rich_au)
    rt = ("\\begin{tabular}{@{}llcccccccc@{}}\\toprule\n"
          "City & & $\\widehat\\ell$ & $\\widehat C$ & $\\widehat r$ & $q_{25}$ & $q_{50}$ & $q_{75}$ & cv & $p_{99}$\\\\ \\midrule\n"
          + "\n".join(rows) + "\n\\bottomrule\\end{tabular}")
    cmd("richTable", rt.replace("\n", " "))
    cmd("richSFMaxZ", f"{rich_sf['gof']['max_abs_z']:.1f}")
    cmd("richAusMaxZ", f"{rich_au['gof']['max_abs_z']:.1f}")
    cmd("richSFHub", f"{rich_sf['held']['hub_share'][0]:.2f}")
    cmd("richSFMaxDeg", f"{rich_sf['held']['max_deg'][0]:.0f}")
    for tag, dd in (("SF", rich_sf), ("Aus", rich_au)):
        cf = dd["counterfactuals"]
        cmd(f"cfHub{tag}", f"{-cf['no_hubs']['delta']:+.3f}")
        cmd(f"cfClo{tag}", f"{-cf['no_closure']['delta']:+.3f}")
        cmd(f"cfLr{tag}", f"{-cf['no_longrange']['delta']:+.3f}")
        cmd(f"cfB{tag}", f"{-cf['b_direct']['delta']:+.3f}")
    zA = rich_au["gof"]["z"]
    cmd("richAusEllZ", f"{zA[0]:+.1f}")
    cmd("richAusPZ", f"{zA[7]:+.1f}")
    others = [abs(zA[i]) for i in (1, 2, 3, 4, 5, 6)]
    cmd("richAusGoodZ", f"{max(others):.1f}")
    tgA = rich_au["target"]["ell"]; mnA = rich_au["gof"]["mean"][0]
    cmd("richAusEllOver", f"{100*(mnA-tgA)/tgA:.0f}")
else:
    for nm in ("richTable richSFMaxZ richAusMaxZ richSFHub richSFMaxDeg "
               "richAusEllZ richAusPZ richAusEllOver "
               "cfHubSF cfCloSF cfLrSF cfBSF cfHubAus cfCloAus cfLrAus "
               "cfBAus richAusGoodZ").split():
        cmd(nm, "---")

if comp and rich_sf:
    cf = rich_sf["counterfactuals"]
    cft = ("\\begin{tabular}{@{}lcc@{}}\\toprule\n"
           "counterfactual (off) & SF Bay & Austin\\\\ \\midrule\n"
           f"propensities $\\sigma$ & {-cf['no_hubs']['delta']:+.3f} & AUSHUB\\\\\n"
           f"closure $\\tau$ & {-cf['no_closure']['delta']:+.3f} & AUSCLO\\\\\n"
           f"long range $w$ & {-cf['no_longrange']['delta']:+.3f} & AUSLR\\\\\n"
           f"dependence set to $\\widehat b_{{\\mathrm{{dir}}}}$ & {-cf['b_direct']['delta']:+.3f} & AUSB\\\\\n"
           "\\bottomrule\\end{tabular}")
    if rich_au:
        cfa = rich_au["counterfactuals"]
        cft = cft.replace("AUSHUB", f"{-cfa['no_hubs']['delta']:+.3f}")
        cft = cft.replace("AUSCLO", f"{-cfa['no_closure']['delta']:+.3f}")
        cft = cft.replace("AUSLR", f"{-cfa['no_longrange']['delta']:+.3f}")
        cft = cft.replace("AUSB", f"{-cfa['b_direct']['delta']:+.3f}")
    obs = comp["observed"]
    def _pp(m, k):
        return comp[m]["pp"][k][0]
    crows = []
    for label, m, aucs in (
            ("geometry ($\\psi{=}0$)", "colas_geo",
             f"{comp['colas_geo']['auc']:.2f}"),
            ("DC-SBM ($K{=}" + str(comp["dcsbm"].get("K", "?")) + "$)",
             "dcsbm", f"{comp['dcsbm']['auc']:.2f}"),
            ("spatial logistic $+$ node effects", "spatial_logistic",
             f"{comp['spatial_logistic']['auc']:.2f}"),
            ("richer class", "richer_class",
             f"{comp['richer_class']['auc']:.2f}/"
             f"{comp['richer_class']['auc_marks_conditional']:.2f}")):
        crows.append(
            f"{label} & {aucs} & {_pp(m,'C'):.3f} & ${_pp(m,'r'):+.3f}$ & "
            f"{_pp(m,'q50'):.1f} & {_pp(m,'deg_cv'):.2f} & {_pp(m,'p99'):.0f}\\\\")
    crows.append(
        f"\\textbf{{observed}} & & \\textbf{{{obs['C']:.3f}}} & "
        f"$\\mathbf{{{obs['r']:+.3f}}}$ & \\textbf{{{obs['q50']:.1f}}} & "
        f"\\textbf{{{obs['deg_cv']:.2f}}} & \\textbf{{{obs['p99']:.0f}}}\\\\")
    ct = ("\\begin{tabular}{@{}lcccccc@{}}\\toprule\n"
          "model & AUC & $C$ & $r$ & $q_{50}$ & cv & $p_{99}$\\\\ \\midrule\n"
          + "\n".join(crows) + "\n\\bottomrule\\end{tabular}")
    cmd("cfCompTable", (cft + "\\\\[7pt]\n" + ct).replace("\n", " "))
    cmd("aucLogit", f"{comp['spatial_logistic']['auc']:.2f}")
    cmd("aucDcsbm", f"{comp['dcsbm']['auc']:.2f}")
    cmd("aucRichCond", f"{comp['richer_class']['auc_marks_conditional']:.2f}")
    cmd("aucGeo", f"{comp['colas_geo']['auc']:.2f}")
    cmd("compLogC", f"{_pp('spatial_logistic','C'):.3f}")
    cmd("compDcsbmC", f"{_pp('dcsbm','C'):.3f}")
else:
    for nm in ("cfCompTable aucLogit aucDcsbm aucRichCond aucGeo compLogC "
               "compDcsbmC").split():
        cmd(nm, "---")

# ---- OpenAlex collaboration network (two-view analysis) ----------------
def _oaj(name):
    p = RUNS + name
    return json.load(open(p)) if os.path.exists(p) else None

oadesc = _oaj("oa_descriptives.json")
oabase = _oaj("oa_base_NLsub.json")
oafree = _oaj("oa_rich_NLsub_free.json")
oafix = _oaj("oa_rich_NLsub_fixdir.json")
oafree2 = _oaj("oa_rich_NLsub_free2.json")

if oadesc:
    rg = oadesc["regions"]
    NL, SUB = rg["Netherlands"], rg["NL-sub"]
    cmd("oaAuthors", "73{,}342"); cmd("oaEdgesAll", "403{,}117")
    cmd("oaNLn", f"{NL['n']:,}".replace(",", "{,}"))
    cmd("oaNLE", f"{NL['edges']:,}".replace(",", "{,}"))
    cmd("oaNLell", f"{NL['ell']:.2f}"); cmd("oaNLC", f"{NL['C']:.3f}")
    cmd("oaNLr", f"{NL['r']:+.3f}")
    cmd("oaNLrs", f"{NL['r_spearman']:+.3f}")
    cmd("oaBdir", f"{NL['b_direct']:+.4f}")
    cmd("oaBdirSE", f"{NL['se_b']:.4f}")
    cmd("oaBdirSEcl", f"{NL['se_b_cluster']:.4f}")
    cmd("oaBdirZcl", f"{NL['b_direct']/NL['se_b_cluster']:.1f}")
    cmd("oaPsiDir", f"{NL['psi_direct']:.4f}")
    cmd("oaNLatoms", str(NL["n_atoms"]))
    cmd("oaNLtop", pct(NL["top_atom_share"], 0))
    cmd("oaGroBdir", f"{rg['Groningen']['b_direct']:+.4f}")
    cmd("oaGroTop", pct(rg["Groningen"]["top_atom_share"], 0))
    cmd("oaSubN", f"{SUB['n']:,}".replace(",", "{,}"))
    cmd("oaSubE", f"{SUB['edges']:,}".replace(",", "{,}"))
    cmd("oaSubEll", f"{SUB['ell']:.2f}"); cmd("oaSubC", f"{SUB['C']:.3f}")
    cmd("oaSubR", f"{SUB['r']:+.3f}")
    cmd("oaSubRs", f"{SUB['r_spearman']:+.3f}")
    cmd("oaSubBdir", f"{SUB['b_direct']:+.4f}")
    cmd("oaSubBdirSEcl", f"{SUB['se_b_cluster']:.4f}")
    cmd("oaSubCfgR", f"{SUB['config_null']['r']:+.3f}")
    cmd("oaSubCfgRsd", f"{SUB['config_null']['r_sd']:.3f}")
    cmd("oaSubCfgC", f"{SUB['config_null']['C']:.4f}")
    q = SUB["edge_km_q"]
    cmd("oaSubQtf", f"{q[0]:.2f}"); cmd("oaSubQf", f"{q[1]:.2f}")
    cmd("oaSubQsf", f"{q[2]:.1f}"); cmd("oaSubQnn", f"{q[3]:.0f}")
    cmd("oaSubFracLocal", pct(SUB["frac_within_1km"], 0))
    cmd("oaSubFracFar", pct(SUB["frac_beyond_6km"], 0))
    cmd("oaSubDbar", f"{SUB['d1_bar']:+.2f}")
    cmd("oaSubDbarZ", f"{SUB['d1_bar']/SUB['d1_se']:.0f}")
    rows = []
    for nm in ("Netherlands", "Randstad", "SouthHolland", "Amsterdam",
               "Utrecht", "Eindhoven", "Groningen"):
        b = rg[nm]
        lab = {"SouthHolland": "South Holland"}.get(nm, nm)
        rows.append(
            f"{lab} & {b['n']:,} & {b['edges']:,} & {b['ell']:.2f} & "
            f"{b['C']:.3f} & ${b['r']:+.3f}$ & ${b['b_direct']:+.4f}$ & "
            f"{b['se_b_cluster']:.4f} & {b['n_atoms']} & "
            f"{pct(b['top_atom_share'],0)}\\\\".replace(",", "{,}")
            .replace("{{,}}", "{,}"))
    cmd("oaRegionTable",
        ("\\begin{tabular}{@{}lrrrrrrrrr@{}}\\toprule\n"
         "region & $n$ & edges & $\\widehat\\ell$ & $\\widehat C$ & "
         "$\\widehat r$ & $\\widehat b_{\\mathrm{dir}}$ & cl.~se & "
         "atoms & top\\\\ \\midrule\n" + "\n".join(rows)
         + "\n\\bottomrule\\end{tabular}").replace("\n", " "))
else:
    for nm in ("oaAuthors oaEdgesAll oaNLn oaNLE oaNLell oaNLC oaNLr oaNLrs "
               "oaBdir oaBdirSE oaBdirSEcl oaBdirZcl oaPsiDir oaNLatoms "
               "oaNLtop oaGroBdir oaGroTop oaSubN oaSubE oaSubEll oaSubC "
               "oaSubR oaSubRs oaSubBdir oaSubBdirSEcl oaSubCfgR oaSubCfgRsd "
               "oaSubCfgC oaSubQtf oaSubQf oaSubQsf oaSubQnn oaSubFracLocal "
               "oaSubFracFar oaSubDbar oaSubDbarZ oaRegionTable").split():
        cmd(nm, "---")

if oabase:
    cmd("oaBaseEll", f"{oabase['model_moments'][0]:.2f}")
    cmd("oaBaseR", f"{oabase['model_moments'][2]:+.3f}")
    cmd("oaBaseGofZ", f"{oabase['gof_z']:.0f}")
    cmd("oaBaseFracOut", pct(oabase["frac_outside"], 1))
    cmd("oaBaseRmax", f"{oabase['R_max_cand']:.2f}")
    cmd("oaBaseDsort", f"{oabase['dsort_signed']:+.3f}")
    cmd("oaBaseDsortFlip", f"{oabase['dsort_flip']:+.3f}")
else:
    for nm in ("oaBaseEll oaBaseR oaBaseGofZ oaBaseFracOut oaBaseRmax "
               "oaBaseDsort oaBaseDsortFlip").split():
        cmd(nm, "---")

def _slope(prof):
    bb = np.array([p[0] for p in prof]); rr = np.array([p[1] for p in prof])
    return float(np.polyfit(bb, rr, 1)[0])

if oafree and oafix:
    zf = oafree["gof"]["z"]; zx = oafix["gof"]["z"]
    cmd("oaFreeB", f"{oafree['theta']['b']:+.2f}")
    cmd("oaFreeMaxZ", f"{oafree['gof']['max_abs_z']:.1f}")
    cmd("oaFreeZell", f"{zf[0]:+.1f}"); cmd("oaFreeZr", f"{zf[2]:+.1f}")
    cmd("oaFreeZp", f"{zf[7]:+.1f}")
    cmd("oaFreeEll", f"{oafree['gof']['mean'][0]:.2f}")
    cmd("oaFixTau", f"{oafix['theta']['tau']:.2f}")
    cmd("oaFixL", f"{oafix['theta']['L']:.0f}")
    cmd("oaFixMaxZ", f"{oafix['gof']['max_abs_z']:.1f}")
    cmd("oaFixZr", f"{zx[2]:+.1f}"); cmd("oaFixZcv", f"{zx[6]:+.1f}")
    cmd("oaFixZqsf", f"{zx[5]:+.1f}")
    cmd("oaFixRhat", f"{oafix['r_hat']:+.3f}")
    cfx = oafix["counterfactuals"]
    cmd("oaFixCfClo", f"{-cfx['no_closure']['delta']:+.3f}")
    cmd("oaFixSortContrib", f"{cfx['b_zero']['delta']:+.4f}")
    cmd("oaFixSignedGap", f"{cfx['b_direct_flip']['delta']:.4f}")
    cmd("oaProfSlopeFix", f"{_slope(oafix['b_profile']):+.3f}")
    cmd("oaProfSlopeFree", f"{_slope(oafree['b_profile']):+.3f}")
    cmd("oaFixSpear", f"{oafix['held']['spearman'][0]:+.3f}")
    ztxt = lambda z: f"{z:+.1f}"
    lrow = []
    lrow.append(f"observed (subsample, $n={oadesc['regions']['NL-sub']['n']:,}$)"
                .replace(",", "{,}")
                + f" & ${oadesc['regions']['NL-sub']['r']:+.3f}$ & "
                f"{oadesc['regions']['NL-sub']['C']:.3f} & & \\\\")
    cn = oadesc["regions"]["NL-sub"]["config_null"]
    lrow.append(f"configuration null & ${cn['r']:+.3f}\\,({cn['r_sd']:.3f})$ & "
                f"{cn['C']:.4f} & degrees & ---\\\\")
    lrow.append(f"hard-disc class & ${oabase['model_moments'][2]:+.3f}$ & "
                f"{oabase['model_moments'][1]:.3f} & "
                f"$(\\ell,C,r)$; misses $\\ell$ & "
                f"$R$, $b$ railed; ${pct(oabase['frac_outside'],0)}\\%$ "
                f"outside $\\widehat R$\\\\")
    lrow.append(f"richer, $b$ free & ${oafree['gof']['mean'][2]:+.3f}$ & "
                f"{oafree['gof']['mean'][1]:.3f} & eight; misses "
                f"$\\ell,r,\\mathrm{{cv}},p_{{99}}$ & $b$ railed "
                f"${oafree['theta']['b']:+.2f}$; $\\max|z|{{=}}"
                f"{oafree['gof']['max_abs_z']:.1f}$\\\\")
    lrow.append(f"richer, $b=\\widehat b_{{\\mathrm{{dir}}}}$ & "
                f"${oafix['gof']['mean'][2]:+.3f}$ & "
                f"{oafix['gof']['mean'][1]:.3f} & eight; misses "
                f"$r,\\mathrm{{cv}}$ & $\\max|z|{{=}}"
                f"{oafix['gof']['max_abs_z']:.1f}$ "
                f"($z_r{{=}}{ztxt(zx[2])}$)\\\\")
    cmd("oaTable",
        ("\\begin{tabular}{@{}lllll@{}}\\toprule\n"
         " & model $r$ & model $C$ & fitted moments & flags\\\\ \\midrule\n"
         + "\n".join(lrow)
         + "\n\\bottomrule\\end{tabular}").replace("\n", " "))
else:
    for nm in ("oaFreeB oaFreeMaxZ oaFreeZell oaFreeZr oaFreeZp oaFreeEll "
               "oaFixTau oaFixL oaFixMaxZ oaFixZr oaFixZcv oaFixZqsf "
               "oaFixRhat oaFixCfClo oaFixSortContrib oaFixSignedGap "
               "oaProfSlopeFix oaProfSlopeFree oaFixSpear oaTable").split():
        cmd(nm, "---")

if oafree2:
    cmd("oaFreeTwoB", f"{oafree2['theta']['b']:+.2f}")
    cmd("oaFreeTwoMaxZ", f"{oafree2['gof']['max_abs_z']:.1f}")
else:
    cmd("oaFreeTwoB", "---"); cmd("oaFreeTwoMaxZ", "---")

oasens = _oaj("oa_sensitivity.json")
if oasens:
    bs = [c["b_direct"] for c in oasens]
    zs = [c["b_direct"] / c["se_b_cluster"] for c in oasens]
    cmd("oaSensLo", f"{min(bs):+.3f}")
    cmd("oaSensHi", f"{max(bs):+.3f}")
    cmd("oaSensZlo", f"{min(zs):.1f}")
    cmd("oaSensZhi", f"{max(zs):.1f}")
else:
    for nm in "oaSensLo oaSensHi oaSensZlo oaSensZhi".split():
        cmd(nm, "---")

with open(OUT, "w") as fh:
    fh.write("% auto-generated from the replication package -- do not edit\n")
    fh.write("\n".join(L) + "\n")
print("rewrote", OUT, f"({len(L)} macros)")

"""Paper figures. Palette: validated slots blue #2a78d6 / orange #eb6834
(+ aqua #1baf7a only with direct labels); ink #0b0b0b, secondary #52514e,
grid #d9d9d9.  Print-safe: series also differ by linestyle/marker."""
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
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
from colas_moments import MomentMap

BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, SEC, GRID = "#0b0b0b", "#52514e", "#d9d9d9"
plt.rcParams.update({
    "font.size": 8.5, "axes.titlesize": 9, "axes.labelsize": 8.5,
    "axes.edgecolor": SEC, "axes.linewidth": 0.8,
    "xtick.color": SEC, "ytick.color": SEC, "text.color": INK,
    "axes.labelcolor": INK, "figure.dpi": 150,
    "axes.spines.top": False, "axes.spines.right": False,
})
RUNS = _p("runs") + _os.sep
FIGS = _p("figures") + _os.sep
R = np.load(RUNS + "sim_results.npz", allow_pickle=True)
mm = MomentMap(N=60, cache=RUNS + "coeffs_N60.pkl")

def style_ax(ax, ygrid=True):
    if ygrid:
        ax.grid(axis="y", color=GRID, lw=0.6, alpha=0.8)
        ax.set_axisbelow(True)

# ================= fig_boundary.pdf (2x2) ==============================
fig, axes = plt.subplots(2, 2, figsize=(6.5, 4.6))

# (a) null distribution of D_n at n=1e4
ax = axes[0, 0]
D = R["null_D_10000"]
mass0 = np.mean(D <= 1e-12)
pos = D[D > 1e-12]
xs = np.linspace(0.02, 8, 300)
bins = np.linspace(0.25, 8, 32)
ax.hist(pos[pos > 0.25], bins=bins, density=False,
        weights=np.full((pos > 0.25).sum(), 1.0 / D.size / (bins[1] - bins[0])),
        color=BLUE, alpha=0.85, edgecolor="white", linewidth=0.4)
frac_small = np.mean((D > 1e-12) & (D <= 0.25)) / 0.25
ax.bar([0.125], [frac_small], width=0.25, color=BLUE, alpha=0.45,
       edgecolor="white", linewidth=0.4)
ax.plot(xs, 0.5 * stats.chi2.pdf(xs, df=1), color=INK, lw=1.6,
        label=r"$\frac{1}{2}\chi^2_1$ density")
ax.bar([-0.55], [mass0], width=0.32, color=ORANGE, edgecolor="white")
ax.annotate(f"mass at 0: {mass0:.3f}\n(theory 0.5)", xy=(-0.5, mass0),
            xytext=(1.6, 0.50), fontsize=7.5, color=SEC,
            arrowprops=dict(arrowstyle="-", color=SEC, lw=0.7))
ax.set_xlim(-0.85, 8); ax.set_ylim(0, 0.62)
ax.annotate("first half-unit\nshown at half tint", xy=(0.35, 0.30),
            fontsize=6.5, color=SEC)
ax.set_xlabel(r"$D_n$"); ax.set_ylabel("density / point mass")
ax.set_title(r"(a) Null distribution of $D_n$,  $n=10^4$", loc="left")
ax.legend(frameon=False, fontsize=7.5)
style_ax(ax)

# (b) QQ of positive part vs chi2_1
ax = axes[0, 1]
q = np.linspace(0.01, 0.995, 120)
emp = np.quantile(pos, q)
theo = stats.chi2.ppf(q, df=1)
ax.plot([0, 9], [0, 9], color=GRID, lw=1.0)
ax.plot(theo, emp, color=BLUE, lw=1.8)
ax.set_xlim(0, 9); ax.set_ylim(0, 9)
ax.set_xlabel(r"$\chi^2_1$ quantile"); ax.set_ylabel(r"empirical quantile of $D_n\,|\,D_n>0$")
ax.set_title("(b) QQ, positive part", loc="left")
style_ax(ax)

# (c) n^{1/4} rate collapse
ax = axes[1, 0]
Ieff = float(R["I_eff_null"])
styles = [(BLUE, "-"), (ORANGE, "--"), (AQUA, ":")]
for (n, (col, ls)) in zip((1000, 10000, 100000), styles):
    psi = R[f"null_psi_{n}"]
    z = n ** 0.25 * np.sqrt(np.maximum(psi, 0.0))
    zs = np.sort(z)
    ax.plot(zs, np.arange(1, zs.size + 1) / zs.size, color=col, ls=ls, lw=1.7,
            label=fr"$n=10^{int(np.log10(n))}$")
t = np.linspace(0, 3.4, 200)
lim = stats.norm.cdf(t ** 2 * np.sqrt(Ieff))
ax.plot(t, lim, color=INK, lw=1.2, alpha=0.9, label="limit law")
ax.set_xlim(0, 3.4); ax.set_ylim(0, 1.02)
ax.set_xlabel(r"$n^{1/4}\,\widehat{b}^{\,+}=n^{1/4}\sqrt{\widehat\psi^{\,+}}$")
ax.set_ylabel("CDF")
ax.set_title(r"(c) $n^{1/4}$ amplitude rate collapse", loc="left")
ax.legend(frameon=False, fontsize=7.5, loc="lower right")
style_ax(ax)

# (d) power curve
ax = axes[1, 1]
psis = R["power_psis"]; pw = R["power"]
ax.axhline(0.05, color=GRID, lw=1.0)
ax.annotate("nominal 5%", xy=(0.06, 0.055), fontsize=7, color=SEC)
ax.plot(psis, pw, color=BLUE, lw=1.8, marker="o", ms=4.5,
        markerfacecolor="white", markeredgewidth=1.4)
for x, y in zip(psis, pw):
    if x in (0.02, 0.04):
        ax.annotate(f"{y:.2f}", xy=(x, y), xytext=(x + 0.002, y - 0.07),
                    fontsize=7, color=SEC)
ax.set_xlim(-0.003, 0.088); ax.set_ylim(0, 1.05)
ax.set_xlabel(r"$\psi$"); ax.set_ylabel("rejection rate")
ax.set_title(r"(d) Power of the 5% boundary test, $n=10^4$", loc="left")
style_ax(ax)

fig.tight_layout(pad=0.8)
fig.savefig(FIGS + "fig_boundary.pdf")
plt.close(fig)
print("fig_boundary.pdf")

# ================= fig_identification.pdf (1x3) ========================
R3 = np.load(RUNS + "sim_r3.npz", allow_pickle=True)
fig, axes = plt.subplots(1, 2, figsize=(6.5, 2.6))

# (a) sigma_min map over (lambda, psi) at eta=9 + certified box
ax = axes[0]
lam = np.linspace(0.6, 5.6, 130); psi = np.linspace(0.0, 0.163, 90)
L, P = np.meshgrid(lam, psi, indexing="ij")
J = mm.targets(L, np.full(L.shape, 9.0), P)["J"]
smin = np.linalg.svd(J, compute_uv=False)[..., -1]
# sequential: one hue, light->dark (blues)
cm = ax.pcolormesh(L, P, smin, cmap="Blues", shading="gouraud",
                   vmin=float(smin.min()) * 0.85, vmax=float(smin.max()),
                   rasterized=True)
cb = fig.colorbar(cm, ax=ax, pad=0.02)
cb.set_label(r"$\sigma_{\min}(DM_3)$ at $\eta=9$", fontsize=7.5)
cb.ax.tick_params(labelsize=7)
ax.add_patch(plt.Rectangle((1, 0), 4, 0.15, fill=False, ec=INK, lw=1.6))
ax.annotate("certified box $B$: global injectivity",
            xy=(1.05, 0.154), fontsize=7.5, color=INK,
            annotation_clip=False)
ax.set_xlabel(r"$\lambda$"); ax.set_ylabel(r"$\psi$")
ax.set_title("(a) Conditioning over $B$", loc="left")

# (b) compensated profiles: r vs psi holding (ell, C) fixed
ax = axes[1]
from scipy.optimize import brentq
anchors = [(2.0, 6.5), (3.0, 9.0), (4.0, 12.0)]
cols = [BLUE, ORANGE, AQUA]
for (lam0, eta0), col in zip(anchors, cols):
    t0 = mm.targets(lam0, eta0, 0.0835)
    ell0, C0 = float(t0["ell"]), float(t0["C"])
    ps_grid = np.linspace(0.0, 0.15, 40)
    rr, ok = [], []
    for ps in ps_grid:
        f = lambda la: float(mm.targets(la, eta0, ps)["C"]) - C0
        try:
            la = brentq(f, 0.4, 6.5, xtol=1e-10)
        except ValueError:
            rr.append(np.nan); continue
        red = mm.reduced(np.array(la), np.array(ps), np.array(ell0))
        rr.append(float(red["rt"]))
    rr = np.array(rr)
    ax.plot(ps_grid, rr, color=col, lw=1.8)
    ax.annotate(fr"$\ell={ell0:.1f},\ C={C0:.2f}$",
                xy=(ps_grid[-1], rr[-1]), xytext=(0.152, rr[-1] - 0.004),
                fontsize=7.5, color=col, ha="left", va="center")
ax.set_xlim(0, 0.21)
ax.set_xticks([0, 0.05, 0.10, 0.15])
ax.set_xlabel(r"$\psi$")
ax.set_ylabel(r"$r$ holding $(\ell, C)$ fixed")
ax.set_title("(b) Compensated profiles", loc="left")
style_ax(ax)

fig.tight_layout(pad=0.8)
fig.savefig(FIGS + "fig_identification.pdf")
plt.close(fig)
print("fig_identification.pdf")

# ================= fig_bias.pdf: O(1/n) validation =====================
fig, ax = plt.subplots(figsize=(3.4, 2.4))
ns = np.array([1000, 10000, 100000])
bias_r = np.array([R[f"moment_bias_{n}"][2] for n in ns])
ax.plot(ns, np.abs(bias_r), color=BLUE, lw=1.8, marker="o", ms=5,
        markerfacecolor="white", markeredgewidth=1.4, label=r"$|\mathrm{bias}(\widehat r_n)|$")
cn = np.abs(bias_r[0]) * ns[0]
ax.plot(ns, cn / ns, color=SEC, lw=1.2, ls="--", label=r"$c/n$")
ax.set_xscale("log"); ax.set_yscale("log")
ax.set_xlabel(r"$n$"); ax.set_ylabel("absolute bias")
ax.set_title(r"Finite-$n$ bias of realized mixing", loc="left")
ax.legend(frameon=False, fontsize=7.5)
style_ax(ax)
fig.tight_layout(pad=0.6)
fig.savefig(FIGS + "fig_bias.pdf")
plt.close(fig)
print("fig_bias.pdf")


# ================= fig_thesis.pdf: the paper's visual thesis ===========
fig, ax = plt.subplots(figsize=(4.6, 3.1))
MA = R3["pairA_m"]; MB = R3["pairB_m"]; shT = R3["pair_shares_true"]
ax.scatter(MA[:, 2], MA[:, 1], s=9, color=BLUE, alpha=0.55, linewidths=0)
ax.scatter(MB[:, 2], MB[:, 1], s=9, color=ORANGE, alpha=0.55, linewidths=0)
ax.annotate(f"sorting share {shT[0]:.2f}", xy=(np.median(MA[:,2]), MA[:,1].mean()+0.011),
            fontsize=9, color=BLUE, ha="center", fontweight="bold")
ax.annotate(f"sorting share ${shT[1]:.2f}$", xy=(np.median(MB[:,2]), MB[:,1].mean()-0.015),
            fontsize=9, color=ORANGE, ha="center", fontweight="bold")
ax.annotate("same mean degree,\nsame assortativity", xy=(0.985, 0.06),
            xycoords="axes fraction", ha="right", fontsize=8, color=SEC)
ax.set_xlabel(r"realized assortativity $\widehat r$  ($n=10^4$)")
ax.set_ylabel(r"realized transitivity $\widehat C$")
style_ax(ax)
fig.tight_layout(pad=0.6)
fig.savefig(FIGS + "fig_thesis.pdf")
plt.close(fig)
print("fig_thesis.pdf")

# ================= fig_annd.pdf: conditional-degree curves =============
import json as _json
SFX = _json.load(open(RUNS + "sf_extras.json")) if __import__("os").path.exists(RUNS + "sf_extras.json") else None
if SFX:
    an = SFX["annd"]
    k = np.array(an["k"]); cnt = np.array(an["count"])
    keep = (k >= 1) & (cnt >= 15)
    fig, ax = plt.subplots(figsize=(4.2, 2.8))
    for key, col, ls, lab in (("obs", BLUE, "-", "Brightkite SF (observed)"),
                              ("fit", INK, "--", "fitted geometry model"),
                              ("cfg", ORANGE, ":", "configuration null")):
        v = np.array(an[key])
        ax.plot(k[keep], v[keep], color=col, ls=ls, lw=1.8, label=lab)
    ax.set_xlabel(r"degree $k$ of one endpoint")
    ax.set_ylabel(r"$\widehat{E}[D_2 \mid D_1=k]$")
    ax.legend(frameon=False, fontsize=7.5)
    style_ax(ax)
    fig.tight_layout(pad=0.6)
    fig.savefig(FIGS + "fig_annd.pdf")
    plt.close(fig)
    print("fig_annd.pdf")

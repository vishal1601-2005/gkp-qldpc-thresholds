"""Regenerate docs/fig*.png. fig1 is computed live (exact Fock-space simulation); fig2/fig3 read data/*.csv."""
import csv, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from gkp_qldpc import fock, sigma2_eff, sigma_from_p, s0sq_from_dB, eta_threshold

# ---- fig1: exact loss channel vs two analytic models ----
etas = (1.0, 0.99, 0.98, 0.95, 0.9, 0.8); fig, ax = plt.subplots(figsize=(6, 4))
for dB, c in ((8, "C0"), (10, "C1"), (12, "C2"), (14, "C3")):
    res = fock.logical_error_vs_loss(dB, etas); xs = [e for e in etas if e < 1]; dr, dp = [], []
    for e in xs:
        p, s0 = res[e]; v = sigma_from_p(p) ** 2
        dr.append(100 * (v / sigma2_eff(s0, e, "rescale") - 1)); dp.append(100 * (v / sigma2_eff(s0, e, "paper") - 1))
    ax.plot(xs, dp, "o--", color=c, alpha=.8, label=f"no-rescale formula, {dB} dB"); ax.plot(xs, dr, "s-", color=c, label=f"corrected, {dB} dB")
ax.axhline(0, color="k", lw=.6); ax.set_xlabel("transmissivity $\\eta$"); ax.set_ylabel("exact-simulation variance vs model (%)")
ax.set_title("Exact Fock-space GKP + loss vs effective-variance models"); ax.invert_xaxis(); ax.legend(fontsize=6, ncol=2); fig.tight_layout(); fig.savefig("docs/fig1_loss_channel_validation.png", dpi=170)

# ---- fig2: block-failure curves, hard vs soft ----
rows = list(csv.DictReader(open("data/layer_a_block_failure_summary.csv"))); fig, axs = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
for ax, dec in zip(axs, ("hard", "soft")):
    for n, c in ((12, "C0"), (20, "C1"), (28, "C2"), (36, "C3"), (40, "C4")):
        pts = [(float(r["sigma"]), float(r["P_block"])) for r in rows if r["decoder"] == dec and int(r["n"]) == n]
        ax.plot(*zip(*pts), "o-", ms=3, color=c, label=f"n={n} (N={n*n+(3*n//4)**2})")
    ax.set_title(f"{dec}-decision BP+OSD"); ax.set_xlabel("GKP shift std $\\sigma$"); ax.grid(alpha=.3)
axs[0].set_ylabel("block failure probability"); axs[0].legend(fontsize=7); fig.tight_layout(); fig.savefig("docs/fig2_code_capacity_curves.png", dpi=170)

# ---- fig3: implied present-day loss under the literal reading of a 'loss indicator' ----
th = list(csv.DictReader(open("data/thresholds_summary.csv"))); s0 = s0sq_from_dB(10)
def implied(sig, conv): return 24.1 * (-10 * np.log10(eta_threshold(sig, s0, conv)))
fig, ax = plt.subplots(figsize=(7, 4)); labels = []
for i, r in enumerate(th):
    lo, hi = float(r["sigma_low"]), float(r["sigma_high"])
    for j, (conv, c) in enumerate((("rescale", "C0"), ("amp", "C3"))):
        a, b = implied(lo, conv), implied(hi, conv)
        ax.plot([a, b], [i + (j - .5) * .28] * 2, "-", lw=6, color=c, solid_capstyle="butt", label=conv if i == 0 else None)
        if abs(b - a) < 1e-9: ax.plot([a], [i + (j - .5) * .28], "D", color=c, ms=7)
    labels.append(f"{r['model']}, {r['decoder']}")
ax.set_yticks(range(len(th))); ax.set_yticklabels(labels, fontsize=8); ax.invert_yaxis()
for v in (3, 10): ax.axvline(v, color="gray", ls=":", lw=1); ax.text(v + .3, len(th) - .55, f"{v} dB", fontsize=7, color="gray")
ax.set_xlabel("implied present-day aggregate loss (dB)\n= 24.1 x threshold loss in dB, at 10 dB squeezing", fontsize=9)
ax.set_title("What a 24.1x loss indicator implies, by modelling choice", fontsize=10); ax.legend(title="loss convention", fontsize=8); fig.tight_layout(); fig.savefig("docs/fig3_implied_loss_sensitivity.png", dpi=170)
for r in th:
    lo, hi = float(r["sigma_low"]), float(r["sigma_high"])
    print(f"{r['model']:26s} {r['decoder']:4s} rescale {implied(lo,'rescale'):5.1f}-{implied(hi,'rescale'):5.1f} dB | amp {implied(lo,'amp'):5.1f}-{implied(hi,'amp'):5.1f} dB")

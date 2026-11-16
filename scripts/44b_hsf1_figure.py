"""Extended Data figure for the HSF1-motif accessibility test (results of scripts/44_hsf1_accessibility.py).

Usage: python -I scripts/44b_hsf1_figure.py
Writes figures/B18_hsf1_accessibility (PDF + PNG): a, HSF1-motif accessibility per donor and glial class in PD, aged
controls and young controls; b, accessibility against the RNA heat-response score of the same donor and class.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.plot.style import apply_style, save_figure, CONDITION_COLORS, INK, MUTED, GRID  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from pdrep import config as cfg  # noqa: E402

B = cfg.RESULTS / "biology"
CLS = ["Oligodendrocytes", "Astrocytes", "Microglia", "OPCs"]


def main():
    apply_style(font__size=7, axes__linewidth=0.6, axes__edgecolor=INK)
    d = pd.read_csv(B / "hsf1_accessibility.tsv", sep="\t")
    t = pd.read_csv(B / "hsf1_accessibility_tests.tsv", sep="\t").set_index("cell_class")
    d["grp"] = np.where(d.disease.eq("PD"), "PD",
                        np.where(d.age_group.astype(str).str.lower().str.contains("young"), "young control", "aged control"))
    cols = {"young control": "#9ad9b6", "aged control": CONDITION_COLORS["Control"], "PD": "#d4a400"}
    fig, axs = plt.subplots(2, 4, figsize=(7.2, 4.0), gridspec_kw={"hspace": 0.65, "wspace": 0.45})
    for j, c in enumerate(CLS):
        x = d[d.cell_class == c]
        ax = axs[0, j]
        for i, g in enumerate(["young control", "aged control", "PD"]):
            v = x[x.grp == g].hsf1_access
            ax.scatter(np.full(len(v), i) + np.random.default_rng(i).uniform(-0.15, 0.15, len(v)), v, s=8,
                       color=cols[g], edgecolor=INK, lw=0.2)
            if len(v):
                ax.plot([i - 0.25, i + 0.25], [v.median()] * 2, color=INK, lw=1.2)
        ax.set_xticks([0, 1, 2], ["young", "aged", "PD"], fontsize=6)
        ax.set_title(f"{c}\nPD vs aged p = {t.loc[c, 'p_pd_vs_ctrl']:.2f}", fontsize=6.5, loc="left")
        if j == 0:
            ax.set_ylabel("HSF1-motif accessibility\n(log ratio to control k-mer)", fontsize=6)
        ax = axs[1, j]
        y = x.dropna(subset=["score"])
        for g in cols:
            v = y[y.grp == g]
            ax.scatter(v.score, v.hsf1_access, s=8, color=cols[g], edgecolor=INK, lw=0.2, label=g)
        ax.set_title(f"ρ = {t.loc[c, 'rho_access_vs_rna_heat']:.2f}, p = {t.loc[c, 'p_rho']:.2f}", fontsize=6.5, loc="left")
        ax.set_xlabel("RNA heat-response score", fontsize=6)
        if j == 0:
            ax.set_ylabel("HSF1-motif accessibility", fontsize=6)
        for a in axs[:, j]:
            for s_ in ["top", "right"]:
                a.spines[s_].set_visible(False)
            a.tick_params(labelsize=5.5, length=2)
    axs[1, 3].legend(frameon=False, fontsize=5.5, loc="lower right")
    fig.text(0.01, 0.98, "a", fontsize=8, fontweight="bold", va="top")
    fig.text(0.01, 0.48, "b", fontsize=8, fontweight="bold", va="top")
    save_figure(fig, ROOT / "figures" / "B18_hsf1_accessibility")
    plt.close(fig)


if __name__ == "__main__":
    main()

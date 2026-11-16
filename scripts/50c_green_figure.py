"""Figure B21: the AD trajectories of Green et al. 2024 against the Lewy-body data.

Usage: python -I scripts/50c_green_figure.py
Inputs from scripts/50a_green_rosmap.py (results/biology/green2024_rosmap.tsv), scripts/50b_green_states_pd.py
(results/biology/green2024_pd[_nocore]/programme_meta.tsv, results/psychad/green2024_specificity.tsv and
green2024_stress_correlation.tsv). Panels: (a) ROSMAP state proportions against AD traits and trajectory (standardised
beta, adjusted for age, sex and PMI); (b) Lewy-body meta-analysis, all signature genes and with heat-shock/core genes
removed (Hedges' g; Lewy stage as Spearman rho); (c) PsychAD DLPFC by diagnosis (control SD, adjusted); (d) PsychAD
donors: correlation of each state score with the cold-inducible score and the glial stress index.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.plot.style import apply_style, save_figure, heatmap, CMAP_DIV, INK  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from pdrep import config as cfg  # noqa: E402

R = cfg.RESULTS
ORDER = ["Mic.12", "Mic.13", "Ast.10", "Oli.7", "Mic.11", "Ast.9", "End.5", "Oli.8", "Ast.5", "Mic.10", "OPC.3", "Mic.6",
         "Mic.15", "Exc.1", "Ast.1", "Mic.14", "Ast.7"]
TRAJ = {"Mic.12": "prAD early", "Mic.13": "prAD early", "Ast.10": "prAD late", "Oli.7": "prAD late",
        "Mic.11": "prAD late, HSP", "Ast.9": "prAD late, HSP", "End.5": "prAD late, HSP", "Oli.8": "ABA, HSP",
        "Ast.5": "ABA", "Mic.10": "ABA", "OPC.3": "ABA", "Mic.6": "ABA", "Mic.15": "ABA", "Exc.1": "ABA, IEG/HSP",
        "Ast.1": "homeostatic", "Mic.14": "homeostatic (IFN)", "Ast.7": "homeostatic (IFN)"}


def state(p):
    return p.split(" ")[0]


def panel(ax, M, F, title, xl, vmax=1.0, fmt="{:.2f}"):
    qm = heatmap(ax, M.values, cmap=CMAP_DIV, vmin=-vmax, vmax=vmax)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            v = M.values[i, j]
            if np.isfinite(v) and F is not None and F.values[i, j] < 0.05:
                ax.text(j, i, "*", ha="center", va="center", fontsize=7, color="white" if abs(v) > 0.6 * vmax else INK)
    ax.set_xticks(range(M.shape[1]), xl, rotation=55, ha="right", fontsize=5.5)
    ax.set_title(title, loc="left", fontsize=7, fontweight="bold")
    return qm


def main():
    apply_style(font__size=6)
    ros = pd.read_csv(R / "biology" / "green2024_rosmap.tsv", sep="\t")
    tr = ["sqrt.amyloid_mf", "sqrt.tangles_mf", "cogng_demog_slope", "prAD", "pmi", "age_death"]
    A = ros.pivot_table(index="state", columns="trait", values="beta").reindex(index=ORDER, columns=tr)
    Af = ros.pivot_table(index="state", columns="trait", values="fdr").reindex(index=ORDER, columns=tr)
    tests = ["PDany_vs_control", "DLB_vs_control", "Lewy stage", "AD_vs_control"]
    pdm = []
    for sub, tag in [("green2024_pd", ""), ("green2024_pd_nocore", " (core removed)")]:
        f = R / "biology" / sub / "programme_meta.tsv"
        if not f.exists():
            continue
        m = pd.read_csv(f, sep="\t")
        m = m[m.test.isin(tests if not tag else ["PDany_vs_control", "Lewy stage"])]
        m["state"] = m.programme.map(state)
        m["col"] = m.test + tag
        pdm.append(m)
    m = pd.concat(pdm)
    m.loc[m.test == "Lewy stage", "est"] = np.tanh(m.loc[m.test == "Lewy stage", "est"])
    cols = [c for c in ["PDany_vs_control", "DLB_vs_control", "AD_vs_control", "Lewy stage",
                        "PDany_vs_control (core removed)", "Lewy stage (core removed)"] if c in set(m.col)]
    B_ = m.pivot_table(index="state", columns="col", values="est").reindex(index=ORDER, columns=cols)
    Bf = m.pivot_table(index="state", columns="col", values="fdr").reindex(index=ORDER, columns=cols)
    ps = pd.read_csv(R / "psychad" / "green2024_specificity.tsv", sep="\t")
    ps["state"] = ps.programme.map(state)
    grp = ["Lewy-body disease", "Lewy-body disease with AD", "Alzheimer's disease", "vascular dementia",
           "frontotemporal dementia", "schizophrenia"]
    n = ps.drop_duplicates("group").set_index("group").n_group
    C = ps.pivot_table(index="state", columns="group", values="effect_sd").reindex(index=ORDER, columns=grp)
    Cf = ps.pivot_table(index="state", columns="group", values="fdr").reindex(index=ORDER, columns=grp)
    co = pd.read_csv(R / "psychad" / "green2024_stress_correlation.tsv", sep="\t")
    co = co[co.donors == "all donors"].assign(state=lambda d: d.programme.map(state))
    D = co.pivot_table(index="state", columns="versus", values="rho").reindex(index=ORDER,
                                                                             columns=["cold_inducible", "stress_index"])
    Df = co.pivot_table(index="state", columns="versus", values="p").reindex(index=ORDER,
                                                                            columns=["cold_inducible", "stress_index"])
    fig, axs = plt.subplots(1, 4, figsize=(7.2, 4.4), gridspec_kw={"width_ratios": [6, len(cols), 6, 2], "wspace": 0.35})
    q = panel(axs[0], A, Af, "a  ROSMAP (Green et al.)",
              ["amyloid", "tangles", "cognitive slope (low = decline)", "prAD fate", "PMI", "age"])
    axs[0].set_yticks(range(len(ORDER)), [f"{s}  ·  {TRAJ[s]}" for s in ORDER], fontsize=5.5)
    lab = {"PDany_vs_control": "PD vs control", "DLB_vs_control": "DLB vs control", "AD_vs_control": "AD vs control",
           "Lewy stage": "Lewy stage (ρ)", "PDany_vs_control (core removed)": "PD vs control, core removed",
           "Lewy stage (core removed)": "Lewy stage, core removed"}
    panel(axs[1], B_, Bf, "b  Lewy-body meta-analysis", [lab[c] for c in cols])
    panel(axs[2], C, Cf, "c  PsychAD DLPFC", [f"{g.replace('Lewy-body disease', 'LBD').replace('frontotemporal', 'FT')} ({n.get(g, 0)})" for g in grp])
    panel(axs[3], D, Df.map(lambda p: p * D.size), "d  PsychAD donors", ["vs CIRBP/RBM3", "vs stress index"])
    for ax in axs[1:]:
        ax.set_yticks([])
    for ax in axs:
        for y in (1.5, 3.5, 7.5, 13.5):
            ax.axhline(y, color="white", lw=1.5)
    cb = fig.colorbar(q, ax=axs, orientation="horizontal", fraction=0.025, pad=0.42, aspect=60)
    cb.set_label("effect (a standardised β; b Hedges' g or ρ; c control SD; d Spearman ρ); * FDR < 0.05",
                 fontsize=5.5)
    fig.subplots_adjust(left=0.17, right=0.99, top=0.93, bottom=0.36)
    save_figure(fig, ROOT / "figures" / "B21_ad_trajectories")
    plt.close(fig)


if __name__ == "__main__":
    main()

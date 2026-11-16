"""Biology figures (main text): composition, consensus transcription, programme map, genetics, evidence map.

Usage: python -I scripts/28_biology_figures.py
Inputs: results/biology/*.tsv(.gz) from scripts 22-27. Writes figures/B2_composition ... B6_evidence (PDF + PNG).
Fig. 1 (study overview) stays figures/F1_landscape; the replication figures F2-F6 become Extended Data.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.plot.style import (apply_style, save_figure, heatmap, CMAP_DIV, ROBUSTNESS_COLORS, ROBUSTNESS_LABELS,  # noqa: E402
                              INK, MUTED, GRID)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch  # noqa: E402

from pdrep import config as cfg  # noqa: E402

B = cfg.RESULTS / "biology"
FIG = ROOT / "figures"
REGION_COLORS = {"all regions": INK, "midbrain / SN": "#2a78d6", "striatum": "#eb6834", "cortex": "#1baf7a"}
UP, DOWN = "#a50026", "#313695"
ROB = ["robust", "supported", "not_replicated", "inconsistent", "contradicted"]


def clean(ax):
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.tick_params(length=2)


def b2():
    m = pd.read_csv(B / "composition_meta.tsv", sep="\t")
    m = m[(m.contrast == "PDany_vs_control")]
    allr = m[(m.denominator == "all nuclei") & (m.scope == "all regions") & (m.k >= 5)].sort_values("est")
    order = list(allr.cell_class)
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.6), gridspec_kw={"width_ratios": [2.2, 1], "wspace": 0.55})
    ax = axs[0]
    off = {"all regions": 0, "midbrain / SN": -0.22, "striatum": 0.22, "cortex": 0.11}
    for scope, col in REGION_COLORS.items():
        d = m[(m.denominator == "all nuclei") & (m.scope == scope) & m.cell_class.isin(order) & (m.k >= 3)]
        for r in d.itertuples():
            y = order.index(r.cell_class) + off[scope]
            big = scope == "all regions"
            ax.errorbar(r.est, y, xerr=1.96 * r.se, fmt="D" if big else "o", color=col, ms=5 if big else 3.2,
                        lw=1.4 if big else 0.8, capsize=0, zorder=3 if big else 2, alpha=1 if big else 0.85)
            if big and r.fdr < 0.05:
                ax.text(r.est + 1.96 * r.se + 0.03, y, f"FDR {r.fdr:.2f}", va="center", fontsize=6, color=INK)
    ax.axvline(0, color=GRID, lw=1, zorder=0)
    ax.set_yticks(range(len(order)), [f"{c} ({int(allr.set_index('cell_class').k[c])})" for c in order])
    ax.set_xlabel("PD vs control: difference in logit share of nuclei (95% CI)")
    ax.set_title("a  Cell-type composition across datasets", loc="left", fontsize=8, fontweight="bold")
    ax.legend(handles=[plt.Line2D([], [], color=c, marker="D" if k == "all regions" else "o", lw=0, ms=4, label=k)
                       for k, c in REGION_COLORS.items()], frameon=False, fontsize=6, loc="lower right")
    clean(ax)
    ax = axs[1]
    d = m[(m.cell_class == "DA neurons") & (m.denominator == "all neurons")]
    e = pd.read_csv(B / "composition_estimates.tsv", sep="\t")
    e = e[(e.contrast == "PDany_vs_control") & (e.cell_class == "DA neurons") & (e.denominator == "all neurons")
          & e.region.isin(["SN", "midbrain"])].sort_values("est")
    for i, r in enumerate(e.itertuples()):
        ax.errorbar(r.est, i, xerr=1.96 * r.se, fmt="o", color="#2a78d6", ms=3.5, lw=0.8)
    pooled = d[d.scope == "midbrain / SN"].iloc[0]
    ax.errorbar(pooled.est, len(e) + 0.6, xerr=1.96 * pooled.se, fmt="D", color=INK, ms=5, lw=1.4)
    ax.set_yticks(list(range(len(e))) + [len(e) + 0.6], list(e.dataset) + [f"pooled (p={pooled.p:.2f})"])
    ax.axvline(0, color=GRID, lw=1, zorder=0)
    ax.set_xlabel("PD vs control: logit share of neurons")
    ax.set_title("b  DA neurons among nigral neurons", loc="left", fontsize=8, fontweight="bold")
    clean(ax)
    save_figure(fig, FIG / "B2_composition")
    plt.close(fig)


def all_regions(d):
    a = d[d.scope == "all regions"].copy()
    a["frac_dir"] = np.where(a.est > 0, a.n_pos / a.k, 1 - a.n_pos / a.k)
    return a


def consensus(d):
    """Consensus genes of a gene meta-analysis (23): all-region scope, FDR < 0.01, |effect| >= 0.2, >= 70% of datasets
    in the pooled direction."""
    a = all_regions(d)
    return a[(a.fdr < 0.01) & (a.est.abs() >= 0.2) & (a.frac_dir >= 0.7)]


def b3():
    d = pd.read_csv(B / "de_meta_PDany_vs_control.tsv.gz", sep="\t")
    a = all_regions(d)
    cons = consensus(d)
    classes = ["Oligodendrocytes", "Astrocytes", "Microglia", "OPCs", "Excitatory neurons", "Inhibitory neurons",
               "Endothelial", "Pericytes"]
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 2.2], width_ratios=[1, 1.4], hspace=0.45, wspace=0.55)
    ax = fig.add_subplot(gs[0, 0])
    n = cons[cons.cell_class.isin(classes)].groupby("cell_class").est.agg(up=lambda e: (e > 0).sum(), down=lambda e: (e < 0).sum())
    n = n.reindex(classes)
    y = np.arange(len(n))
    ax.barh(y, n.up, color=UP, height=0.38, label="up in PD")
    ax.barh(y + 0.4, n.down, color=DOWN, height=0.38, label="down in PD")
    ax.set_yticks(y + 0.2, n.index)
    ax.invert_yaxis()
    ax.set_xlabel("consensus genes")
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    ax.set_title("a  Consensus PD genes per cell type", loc="left", fontsize=8, fontweight="bold")
    clean(ax)
    # shared genes: consensus in >= 4 of the 8 classes
    c8 = cons[cons.cell_class.isin(classes)]
    shared = c8.groupby("gene").agg(n=("cell_class", "nunique"), z=("z", "mean"))
    top_up = shared[(shared.n >= 4) & (shared.z > 0)].sort_values(["n", "z"], ascending=False).head(22).index
    top_dn = shared[(shared.n >= 4) & (shared.z < 0)].sort_values(["n", "z"], ascending=[False, True]).head(8).index
    genes = list(top_up) + list(top_dn)
    M = a[a.cell_class.isin(classes) & a.gene.isin(genes)].pivot_table(index="gene", columns="cell_class", values="est")
    M = M.reindex(index=genes, columns=classes)
    ax = fig.add_subplot(gs[:, 1])
    qm = heatmap(ax, M.values, cmap=CMAP_DIV, vmin=-1, vmax=1)
    ax.set_xticks(range(len(classes)), classes, rotation=90)
    ax.set_yticks(range(len(genes)), genes, fontsize=6)
    sig = c8.set_index(["gene", "cell_class"]).index
    for i, g in enumerate(genes):
        for j, c in enumerate(classes):
            if (g, c) in sig:
                ax.text(j, i, "•", ha="center", va="center", fontsize=6, color="white")
    cb = fig.colorbar(qm, ax=ax, fraction=0.04, pad=0.03)
    cb.set_label("pooled log2 fold change", fontsize=6)
    cb.ax.tick_params(labelsize=6)
    ax.set_title("c  Genes changed in PD across cell types", loc="left", fontsize=8, fontweight="bold")
    ax = fig.add_subplot(gs[1, 0])
    g = pd.read_csv(B / "gsea_PDany_vs_control.tsv.gz", sep="\t") if (B / "gsea_PDany_vs_control.tsv.gz").exists() else None
    if g is not None:
        g = g[(g.scope == "all regions") & g.cell_class.isin(classes)]
        h = g[g.library == "Hallmark"]
        keep = list(h[h.fdr < 0.05].groupby("term").size().sort_values(ascending=False).head(13).index)
        # the strongest Reactome signals (heat-shock factor and integrated stress response), shown on top
        rea = ["REACTOME_HSF1_ACTIVATION", "REACTOME_ATTENUATION_PHASE", "REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS",
               "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY"]
        sub = g[g.term.isin(keep + rea)]
        P = sub.pivot_table(index="term", columns="cell_class", values="nes").reindex(columns=classes)
        F = sub.pivot_table(index="term", columns="cell_class", values="fdr").reindex(columns=classes)
        hm = P.loc[keep]
        order = rea + list(hm.mean(1).sort_values(ascending=False).index)
        P, F = P.loc[order], F.loc[order]
        names = {"REACTOME_HSF1_ACTIVATION": "HSF1 activation (R)", "REACTOME_ATTENUATION_PHASE": "HSF1 attenuation phase (R)",
                 "REACTOME_CELLULAR_RESPONSE_TO_HEAT_STRESS": "response to heat stress (R)",
                 "REACTOME_RESPONSE_OF_EIF2AK4_GCN2_TO_AMINO_ACID_DEFICIENCY": "GCN2 / amino-acid deficiency (R)"}
        P.index = [names.get(t, t.replace("HALLMARK_", "").replace("_", " ").lower()) for t in P.index]
        qm = heatmap(ax, P.values, cmap=CMAP_DIV, vmin=-2.5, vmax=2.5)
        for i in range(P.shape[0]):
            for j in range(P.shape[1]):
                if F.values[i, j] < 0.05:
                    ax.text(j, i, "•", ha="center", va="center", fontsize=6, color="white")
        ax.axhline(len(rea) - 0.5, color="white", lw=1.5)
        ax.set_xticks(range(len(classes)), classes, rotation=90)
        ax.set_yticks(range(len(P)), P.index, fontsize=6)
        cb = fig.colorbar(qm, ax=ax, fraction=0.05, pad=0.03)
        cb.set_label("NES", fontsize=6)
        cb.ax.tick_params(labelsize=6)
    else:  # placeholder until the GSEA run (scripts/24_pathways.py) has been merged
        ax.text(0.5, 0.5, "Pathway enrichment\n(preranked GSEA)\n\nresults pending", ha="center", va="center",
                fontsize=8, color=MUTED, transform=ax.transAxes)
        ax.set_xticks([])
        ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color(GRID)
    ax.set_title("b  Pathways (preranked GSEA)", loc="left", fontsize=8, fontweight="bold")
    save_figure(fig, FIG / "B3_consensus_transcription")
    plt.close(fig)
    cons.to_csv(B / "consensus_genes_PDany.tsv", sep="\t", index=False)


def b4():
    m = pd.read_csv(B / "programme_meta.tsv", sep="\t")
    tests = ["PD_vs_control", "PDD_vs_control", "DLB_vs_control", "LBearly_vs_control", "LBmid_vs_control",
             "LBlate_vs_control", "Lewy stage", "AD_vs_control", "post-mortem interval", "age"]
    labels = ["PD", "PDD", "DLB", "Lewy early\n(Braak 1–2)", "Lewy mid\n(3–4)", "Lewy late\n(5–6)", "Lewy stage\n(trend)",
              "AD", "post-mortem\ninterval", "age"]
    general = [("heat response", "Oligodendrocytes"), ("heat response", "Astrocytes"), ("heat response", "Microglia"),
               ("chaperone folding", "Oligodendrocytes"), ("unfolded protein response", "Oligodendrocytes"),
               ("hypoxia", "Microglia"), ("cold-inducible (CIRBP, RBM3)", "ALL")]
    # published glial states, grouped by cell type; studies interleaved, stress-like states first
    states = [("UPR astrocytes, Astrocytes3 (Martirosyan)", "Astrocytes"), ("UPR astrocytes (Barba-Reyes)", "Astrocytes"),
              ("GFAP/MT astrocytes, Astrocytes4 (Martirosyan)", "Astrocytes"), ("activated astrocytes, Ast-1 (Xu)", "Astrocytes"),
              ("C3/CD44 astrocytes, Astrocytes1 (Martirosyan)", "Astrocytes"), ("AR2 astrocytes (Barba-Reyes)", "Astrocytes"),
              ("reactive astrocytes (Barba-Reyes)", "Astrocytes"),
              ("UPR microglia (Barba-Reyes)", "Microglia"), ("HSP90AA1-high microglia (Smajić)", "Microglia"),
              ("GPNMB-high microglia (Smajić)", "Microglia"), ("activated microglia, Micr-1 (Xu)", "Microglia"),
              ("activated myeloid cells (Ma)", "Microglia"), ("DAM microglia (Barba-Reyes)", "Microglia"),
              ("reactive microglia (Barba-Reyes)", "Microglia"),
              ("heat-shock oligodendroglia (Barba-Reyes)", "Oligodendrocytes"),
              ("disease-associated oligodendrocytes (Adams)", "Oligodendrocytes"),
              ("stress oligodendrocytes, Oligos5 (Martirosyan)", "Oligodendrocytes")]
    other = [("inflammatory response", "Microglia"), ("interferon gamma", "Microglia"), ("oxidative phosphorylation", "ALL")]
    piv = m.set_index(["programme", "cell_class", "test"])
    rows = [r for r in general + states + other if (r[0], r[1], "PD_vs_control") in piv.index]
    Z = np.array([[piv.z.get((p, c, t), np.nan) for t in tests] for p, c in rows])
    F = np.array([[piv.fdr.get((p, c, t), np.nan) for t in tests] for p, c in rows])
    fig, ax = plt.subplots(figsize=(7.2, 7.4))
    qm = heatmap(ax, np.clip(Z, -6, 6), cmap=CMAP_DIV, vmin=-6, vmax=6)
    for i in range(len(rows)):
        for j in range(len(tests)):
            if np.isfinite(F[i, j]) and F[i, j] < 0.05:
                ax.text(j, i, "*", ha="center", va="center", fontsize=8, color="white" if abs(Z[i, j]) > 3.5 else INK)
    ax.set_xticks(range(len(tests)), labels, fontsize=6)
    ax.set_yticks(range(len(rows)), [f"{p} · {c.lower() if c != 'ALL' else 'all cells'}" for p, c in rows], fontsize=6.3)
    ax.axvline(6.5, color="white", lw=2)
    ax.axvline(7.5, color="white", lw=2)
    ng = sum(r in rows for r in general)
    ns = sum(r in rows for r in states)
    for y in [ng - 0.5, ng + ns - 0.5]:
        ax.axhline(y, color="white", lw=2.5)
    for k in range(ng, ng + ns - 1):  # thin lines between cell types within the published states
        if rows[k][1] != rows[k + 1][1]:
            ax.axhline(k + 0.5, color="white", lw=1)
    cb = fig.colorbar(qm, ax=ax, fraction=0.03, pad=0.02)
    cb.set_label("meta-analysis z (vs control, or correlation)", fontsize=6)
    cb.ax.tick_params(labelsize=6)
    ax.set_title("Disease programmes and published glial states: what they track across datasets", loc="left",
                 fontsize=8, fontweight="bold")
    save_figure(fig, FIG / "B4_programme_map")
    plt.close(fig)


def b5():
    g = pd.read_csv(B / "genetics_meta.tsv", sep="\t")
    panels = [("specificity", "a  Common-variant risk (MAGMA)\nin cell-type-specific genes"),
              ("PD-up genes (PDany_vs_control)", "b  Common-variant risk in genes\nup-regulated in PD"),
              ("familial PD genes", "c  Familial PD genes in\ncell-type-specific genes")]
    fig, axs = plt.subplots(1, 3, figsize=(7.2, 3.3), gridspec_kw={"wspace": 0.95})
    for ax, (st, title) in zip(axs, panels):
        d = g[(g.set_type == st) & (g.k >= 3) & (g.cell_class != "ALL")].copy()
        if st == "familial PD genes":  # discrete overlap tests: datasets with p < 0.05 out of those tested
            d["nsig"] = (d.share_p05 * d.k).round()
            d = d.sort_values(["nsig", "k"])
            ax.barh(range(len(d)), d.k, color="#e6e6e6", height=0.65, label="datasets tested")
            ax.barh(range(len(d)), d.nsig, color=[UP if p < 0.05 else "#9e9e9e" for p in d.p], height=0.65,
                    label="datasets with p < 0.05")
            for i, r in enumerate(d.itertuples()):
                if r.p < 0.05:
                    ax.text(r.k + 0.3, i, f"binomial p = {r.p:.3f}", va="center", fontsize=6, color=INK)
            ax.set_xlabel("datasets", fontsize=6.5)
            ax.legend(frameon=False, fontsize=5.5, loc="center right")
        else:
            d["z"] = d.z.clip(-4, 6)
            d = d.sort_values("z")
            ax.barh(range(len(d)), d.z, color=[UP if p < 0.05 else "#bdbdbd" for p in d.p], height=0.65)
            ax.axvline(1.645, color=MUTED, lw=0.8, ls=":")
            ax.axvline(0, color=GRID, lw=1)
            ax.set_xlabel("combined z across datasets", fontsize=6.5)
        ax.set_yticks(range(len(d)), [f"{c} ({k})" for c, k in zip(d.cell_class, d.k)], fontsize=6)
        ax.set_title(title, loc="left", fontsize=7.5, fontweight="bold")
        clean(ax)
    save_figure(fig, FIG / "B5_genetics")
    plt.close(fig)


def b6():
    em = pd.read_csv(B / "evidence_map.tsv", sep="\t", index_col=0)
    em = em.sort_values("share_robust_or_supported")
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    left = np.zeros(len(em))
    for r in ROB:
        v = em.get(r, pd.Series(0, index=em.index)).fillna(0) / em.tested
        ax.barh(range(len(em)), v, left=left, color=ROBUSTNESS_COLORS[r], edgecolor="white", linewidth=1, height=0.7,
                label=ROBUSTNESS_LABELS[r])
        left += v.values
    for i, (t, r) in enumerate(em.iterrows()):
        ax.text(1.01, i, f"{int(r.tested)} claims · {int(r.papers)} papers", va="center", fontsize=6, color=MUTED)
    ax.set_yticks(range(len(em)), em.index)
    ax.set_xlim(0, 1.32)
    ax.set_xticks([0, 0.5, 1], ["0%", "50%", "100%"])
    ax.set_xlabel("share of the theme's claims tested in ≥ 2 other datasets")
    ax.legend(frameon=False, fontsize=6, loc="upper center", bbox_to_anchor=(0.4, -0.18), ncol=5)
    ax.set_title("How published claims hold up, by biological theme", loc="left", fontsize=8, fontweight="bold")
    clean(ax)
    save_figure(fig, FIG / "B6_evidence_map")
    plt.close(fig)


def main():
    apply_style(font__size=7, axes__titlesize=8, axes__labelsize=7, xtick__labelsize=6.5, ytick__labelsize=6.5,
                axes__linewidth=0.6, axes__edgecolor=INK)
    for f in (b2, b3, b4, b5, b6):
        f()
    print("biology figures written")


if __name__ == "__main__":
    main()

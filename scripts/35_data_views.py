"""The data behind the pooled results: donor-level expression of core genes and cell shares in every dataset.

Usage: python -I scripts/35_data_views.py
Reads the pseudobulk counts (results/de_cache/<ds>/input/<region>__<class>.csv.gz, units.tsv) and
results/biology/donor_composition.tsv.gz. Writes results/biology/donor_genes.tsv.gz (log2 CPM of core genes per donor
unit and class) and figures/B10_data_views (PDF + PNG):
  a, b  HSPA1A and CIRBP in oligodendrocytes per donor unit, every dataset x region (relative to the control median)
  c     microglia share of nuclei per donor unit
  d     case-control difference in mean log2 CPM of core genes in oligodendrocytes, per dataset x region
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.plot.style import apply_style, save_figure, heatmap, CMAP_DIV, CONDITION_COLORS, INK, MUTED, GRID, REGION_GROUP  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

from pdrep import config as cfg  # noqa: E402

BASE = cfg.RESULTS
DE, B, FIG = BASE / "de_cache", BASE / "biology", ROOT / "figures"
CASES = {"PD", "PDD", "DLB", "LBD"}
GENES = ["HSPA1A", "HSPA1B", "HSPH1", "DNAJB1", "HSP90AA1", "BAG3", "CHORDC1", "SLC38A2", "P4HA1", "BNIP3L", "CIRBP", "RBM3"]
CLASSES = {"Oligodendrocytes", "Astrocytes", "Microglia", "ALL"}
GROUP_ORDER = ["midbrain / SN", "striatum", "cortex", "other"]
CASE_C, CTRL_C = "#d4a400", CONDITION_COLORS["Control"]
LABEL = {"adams2024": "Adams 2024", "barbareyes2026": "Barba-Reyes 2026", "chia2025": "Chia 2025",
         "dehestani2024": "Dehestani 2024", "feleke2021": "Feleke 2021", "glock2026": "Glock 2026", "jin2024": "Jin 2024",
         "kamath2022": "Kamath 2022", "lee2023": "Lee 2023", "lee2026": "Lee 2026", "ma2025": "Ma 2025",
         "martirosyan2024": "Martirosyan 2024", "menon2026": "Menon 2026", "prashant2024": "Prashant 2024",
         "smajic2022": "Smajić 2022", "volpato2025": "Volpato 2025", "wang2024": "Wang 2024", "xu2023": "Xu 2023",
         "zhu2024": "Zhu 2024"}


def collect():
    rows = []
    for ds_dir in sorted(DE.iterdir()):
        uf = ds_dir / "input" / "units.tsv"
        if not uf.exists():
            continue
        units = pd.read_csv(uf, sep="\t", dtype=str).set_index("unit")
        for f in sorted((ds_dir / "input").glob("*__*.csv.gz")):
            region, cls = f.name[: -len(".csv.gz")].split("__", 1)
            cls = cls.replace("_", " ")
            if cls not in CLASSES:
                continue
            x = pd.read_csv(f, index_col=0)
            x.index = x.index.str.upper()
            lc = np.log2(x.div(x.sum(0), axis=1) * 1e6 + 1)
            g = [k for k in GENES if k in lc.index]
            d = lc.loc[g].T
            d["dataset"], d["region"], d["cell_class"] = ds_dir.name, region, cls
            d["disease"] = units.loc[d.index, "disease"].values
            d["donor"] = units.loc[d.index, "donor"].values
            rows.append(d.rename_axis("unit").reset_index())
    out = pd.concat(rows, ignore_index=True)
    out.to_csv(B / "donor_genes.tsv.gz", sep="\t", index=False)
    return out


def order_units(d):
    u = d[["dataset", "region"]].drop_duplicates().copy()
    u["grp"] = u.region.map(REGION_GROUP).fillna("other")
    u["g"] = u.grp.map({k: i for i, k in enumerate(GROUP_ORDER)}).fillna(3)
    return u.sort_values(["g", "dataset", "region"]).reset_index(drop=True)


def strip(ax, d, col, units, xlabel, center=True):
    for i, r in enumerate(units.itertuples()):
        x = d[(d.dataset == r.dataset) & (d.region == r.region)]
        ref = x[x.disease.eq("control")][col]
        if len(ref) < 2:
            continue
        c0 = ref.median() if center else 0
        for grp, off, colr in [("control", -0.17, CTRL_C), ("case", 0.17, CASE_C)]:
            v = (x[x.disease.eq("control")] if grp == "control" else x[x.disease.isin(CASES)])[col] - c0
            jit = np.random.default_rng(i).uniform(-0.09, 0.09, len(v))
            ax.scatter(v, i + off + jit, s=3, color=colr, alpha=0.75, lw=0, rasterized=True)
            if len(v):
                ax.plot([v.median()] * 2, [i + off - 0.15, i + off + 0.15], color=INK if grp == "case" else "#118a50", lw=1.2)
    ax.set_yticks(range(len(units)), [f"{LABEL.get(r.dataset, r.dataset)} · {r.region}" for r in units.itertuples()],
                  fontsize=5.3)
    ax.set_ylim(len(units) - 0.4, -0.6)
    ax.axvline(0, color=GRID, lw=0.8, zorder=0)
    ax.set_xlabel(xlabel)
    for s in ["top", "right"]:
        ax.spines[s].set_visible(False)
    ax.tick_params(length=2)
    prev = None
    for i, r in enumerate(units.itertuples()):  # region-group separators
        if prev is not None and r.grp != prev:
            ax.axhline(i - 0.5, color=MUTED, lw=0.4, ls=":")
        prev = r.grp


def main():
    apply_style(font__size=7, axes__titlesize=8, axes__labelsize=7, xtick__labelsize=6.5, ytick__labelsize=6.5,
                axes__linewidth=0.6, axes__edgecolor=INK)
    g = collect()
    ol = g[g.cell_class.eq("Oligodendrocytes")]
    units = order_units(ol)
    fig = plt.figure(figsize=(7.2, 9.6))
    gs = fig.add_gridspec(2, 3, height_ratios=[1.35, 1], hspace=0.28, wspace=0.12, width_ratios=[1, 1, 1])
    ax = fig.add_subplot(gs[0, 0])
    strip(ax, ol, "HSPA1A", units, "HSPA1A, oligodendrocytes\n(log2 CPM − control median)")
    ax.set_title("a  HSPA1A in every dataset", loc="left", fontsize=8, fontweight="bold")
    ax = fig.add_subplot(gs[0, 1])
    strip(ax, ol, "CIRBP", units, "CIRBP, oligodendrocytes\n(log2 CPM − control median)")
    ax.set_yticklabels([])
    ax.set_title("b  CIRBP in every dataset", loc="left", fontsize=8, fontweight="bold")
    ax = fig.add_subplot(gs[0, 2])
    comp = pd.read_csv(B / "donor_composition.tsv.gz", sep="\t")
    comp = comp[comp.cell_class.eq("Microglia")].copy()
    comp[["donor", "region"]] = comp.unit.str.split("|", n=1, expand=True)
    dis = g.drop_duplicates(["dataset", "donor"]).set_index(["dataset", "donor"]).disease
    comp["disease"] = [dis.get((a_, b_), np.nan) for a_, b_ in zip(comp.dataset, comp.donor)]
    comp["pct"] = 100 * comp.share
    strip(ax, comp, "pct", units, "microglia\n(% of nuclei)", center=False)
    ax.set_yticklabels([])
    ax.set_xlim(left=0)
    ax.set_title("c  Microglia share", loc="left", fontsize=8, fontweight="bold")
    ax.scatter([], [], s=8, color=CTRL_C, label="controls")
    ax.scatter([], [], s=8, color=CASE_C, label="Lewy-body cases")
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    # d: case-control differences of core genes per dataset x region
    ax = fig.add_subplot(gs[1, :])
    rows = []
    for r in units.itertuples():
        x = ol[(ol.dataset == r.dataset) & (ol.region == r.region)]
        a_, c_ = x[x.disease.isin(CASES)], x[x.disease.eq("control")]
        if len(a_) >= 2 and len(c_) >= 2:
            rows.append(pd.Series(a_[GENES].mean() - c_[GENES].mean(),
                                  name=f"{LABEL.get(r.dataset, r.dataset)} · {r.region} ({len(a_)}/{len(c_)})"))
    M = pd.DataFrame(rows)
    qm = heatmap(ax, M.T.values, cmap=CMAP_DIV, vmin=-1.2, vmax=1.2)
    ax.set_xticks(range(len(M)), M.index, rotation=90, fontsize=5.3)
    ax.set_yticks(range(len(GENES)), GENES, fontsize=6.5)
    cb = fig.colorbar(qm, ax=ax, fraction=0.02, pad=0.01)
    cb.set_label("cases − controls\n(mean log2 CPM)", fontsize=6)
    cb.ax.tick_params(labelsize=6)
    up = (M[GENES[:10]] > 0).mean().mean()
    dn = (M[["CIRBP", "RBM3"]] < 0).mean().mean()
    ax.set_title("d  Core genes in oligodendrocytes, every dataset × region (cases/controls)", loc="left", fontsize=8,
                 fontweight="bold")
    ax.text(0, -0.09, f"{up:.0%} of stress-gene and {dn:.0%} of cold-gene differences point in the pooled direction",
            transform=ax.transAxes, fontsize=6, color=MUTED, va="bottom")
    save_figure(fig, FIG / "B10_data_views")
    plt.close(fig)
    print(len(units), "dataset x region units;", round(up, 2), round(dn, 2))


if __name__ == "__main__":
    import warnings
    warnings.simplefilter("ignore")
    main()

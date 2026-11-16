"""Brain map of the stress response: case-control effects pooled per region and drawn on a schematic brain.

Usage: python -I scripts/39_region_map.py
Input: results/biology/donor_regions_effects.tsv (Hedges' g of cases vs controls per dataset x region x class x
programme; scripts/30_donor_analyses.py). Regions are merged into 12 anatomical sites (SN with midbrain, PFC with
DLPFC, MTG with temporal cortex), pooled per site across datasets by DerSimonian-Laird and drawn at schematic positions
on a midsagittal outline (deep nuclei and lateral cortex projected onto it). Circle colour gives the pooled g, size the
number of donors. Writes results/biology/region_map.tsv and figures/B12_region_map (PDF + PNG).
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.meta import dersimonian_laird  # noqa: E402
from pdrep.plot.style import apply_style, save_figure, CMAP_DIV, INK, MUTED, GRID  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Ellipse, Polygon  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402

from pdrep import config as cfg  # noqa: E402

B = cfg.RESULTS / "biology"
SITE = {"SN": "substantia nigra", "midbrain": "substantia nigra", "PFC": "prefrontal cortex", "DLPFC": "prefrontal cortex",
        "MTG": "temporal cortex", "temporal": "temporal cortex", "cingulate": "cingulate cortex", "amygdala": "amygdala",
        "putamen": "putamen", "caudate": "caudate", "GPi": "globus pallidus", "motor": "motor cortex",
        "visual": "visual cortex", "cerebellum": "cerebellum", "DMNX": "vagal nucleus (DMNX)"}
# schematic midsagittal positions (front to the left); lateral and deep structures projected
POS = {"prefrontal cortex": (1.6, 6.1), "motor cortex": (4.8, 7.55), "cingulate cortex": (4.4, 5.8),
       "visual cortex": (8.6, 5.0), "temporal cortex": (3.9, 2.75), "amygdala": (2.8, 3.45), "caudate": (3.3, 5.0),
       "putamen": (3.85, 4.3), "globus pallidus": (4.45, 3.95), "substantia nigra": (5.35, 3.25),
       "vagal nucleus (DMNX)": (6.0, 1.15), "cerebellum": (7.7, 2.0)}
# label offset (dx, dy) and alignment per site, to keep labels apart
LAB = {"prefrontal cortex": (-0.1, 0.45, "center", "bottom"), "motor cortex": (0, 0.45, "center", "bottom"),
       "cingulate cortex": (0.4, 0.15, "left", "center"), "visual cortex": (0, 0.4, "center", "bottom"),
       "temporal cortex": (0, -0.4, "center", "top"), "amygdala": (-0.35, 0, "right", "center"),
       "caudate": (-0.35, 0.1, "right", "center"), "putamen": (-0.35, 0.05, "right", "center"),
       "globus pallidus": (0.35, 0.3, "left", "bottom"), "substantia nigra": (0.4, 0, "left", "center"),
       "vagal nucleus (DMNX)": (0.4, 0, "left", "center"), "cerebellum": (0.35, 0.3, "left", "bottom")}
BRAAK = {"vagal nucleus (DMNX)": 1, "substantia nigra": 3, "amygdala": 4, "cingulate cortex": 4, "temporal cortex": 4,
         "prefrontal cortex": 5, "motor cortex": 6, "visual cortex": 6}
PANELS = [("heat response", "ALL", "a  Heat response, all cells"), ("heat response", "Oligodendrocytes", "b  Heat response, oligodendrocytes"),
          ("chaperone folding", "Microglia", "c  Chaperone folding, microglia"),
          ("cold-inducible (CIRBP, RBM3)", "ALL", "d  CIRBP / RBM3, all cells")]


def pool_sites(e):
    """Per dataset x region effects (30) -> per programme x cell class x brain site, DerSimonian-Laird."""
    e = e.copy()
    e["site"] = e.region.map(SITE)
    rows = []
    for (p, c, s), d in e.dropna(subset=["site"]).groupby(["programme", "cell_class", "site"]):
        r = dersimonian_laird(d.hedges_g, d.se)
        rows.append({"programme": p, "cell_class": c, "site": s, **r, "donors": int((d.n_cases + d.n_ctrl).sum()),
                     "datasets": ";".join(sorted(d.dataset.unique()))})
    return pd.DataFrame(rows)


def pooled():
    out = pool_sites(pd.read_csv(B / "donor_regions_effects.tsv", sep="\t"))
    out.to_csv(B / "region_map.tsv", sep="\t", index=False)
    return out


def brain(ax):
    t = np.linspace(0, np.pi, 80)
    top = np.c_[5 - 4.6 * np.cos(t), 4.3 + 3.9 * np.sin(t)]  # cerebrum dome
    base = np.array([[9.6, 4.3], [9.2, 3.3], [7.8, 3.05], [6.4, 3.3], [5.6, 2.2], [3.4, 2.1], [1.3, 3.0], [0.4, 4.3]])
    ax.add_patch(Polygon(np.r_[top, base], closed=True, fc="#f3f4f2", ec=MUTED, lw=0.8, zorder=0))
    ax.add_patch(Ellipse((7.7, 2.0), 3.0, 1.9, fc="#ecedea", ec=MUTED, lw=0.8, zorder=0))  # cerebellum
    ax.add_patch(Polygon([[5.0, 3.2], [6.0, 3.2], [6.5, 0.2], [5.6, 0.2]], closed=True, fc="#ecedea", ec=MUTED, lw=0.8,
                         zorder=0))  # brainstem
    ax.plot([2.4, 3.6, 5.6, 7.2], [5.05, 5.35, 5.35, 4.9], color=GRID, lw=2.2, zorder=0, solid_capstyle="round")  # callosum
    ax.set_xlim(-0.9, 10.6)
    ax.set_ylim(-0.4, 8.6)
    ax.set_aspect("equal")
    ax.axis("off")


def main():
    apply_style(font__size=7, axes__titlesize=8)
    m = pooled()
    norm = Normalize(-0.8, 0.8)
    cmap = plt.get_cmap(CMAP_DIV)
    fig, axs = plt.subplots(2, 2, figsize=(7.2, 5.6))
    for ax, (p, c, title) in zip(axs.ravel(), PANELS):
        brain(ax)
        d = m[(m.programme == p) & (m.cell_class == c)]
        for r in d.itertuples():
            x, y = POS[r.site]
            ax.scatter(x, y, s=12 + 7 * np.sqrt(r.donors), color=cmap(norm(r.est)), edgecolor=INK if r.p < 0.05 else MUTED,
                       lw=1.1 if r.p < 0.05 else 0.5, zorder=3)
            lab = r.site.replace(" cortex", " ctx").replace("vagal nucleus (DMNX)", "DMNX").replace("substantia nigra", "SN").replace("globus pallidus", "GPi")
            stage = f" (B{BRAAK[r.site]})" if r.site in BRAAK else ""
            dx, dy, ha, va = LAB[r.site]
            ax.text(x + dx, y + dy, f"{lab}{stage} {r.est:+.2f}", ha=ha, va=va, fontsize=5.2, color=INK, zorder=4)
        ax.set_title(title, loc="left", fontsize=8, fontweight="bold")
    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cb = fig.colorbar(sm, ax=axs, fraction=0.025, pad=0.02)
    cb.set_label("cases vs controls, pooled Hedges' g", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    fig.text(0.02, 0.01, "Circle size: donors; black outline: p < 0.05. B1–B6: Braak Lewy stage at which a region is "
             "typically reached (striatum, pallidum and cerebellum are not part of the Braak sequence).",
             fontsize=5.8, color=MUTED)
    save_figure(fig, ROOT / "figures" / "B12_region_map")
    plt.close(fig)
    print(m[m.cell_class.eq("ALL") & m.programme.eq("heat response")][["site", "est", "p", "k", "donors"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()

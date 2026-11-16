"""Figure style for the PD claim-replication manuscript.

Copied from pd_striatum/revision/scripts/style.py (2026-10-23) so figures match the GLIAD manuscript; additions:
OUTCOME_COLORS / ROBUSTNESS_COLORS (diverging, sampled to pass the dataviz CVD and normal-vision checks on white:
navy / mid-blue / neutral grey / orange / red) and CLAIM_TYPE_COLORS. Original notes follow.


Matches the manuscript figures:
* font: Arial (requested 2026-10-14). The TTF files live in revision/env/fonts/ (Arial is not installed on the
  cluster) and are registered with matplotlib by apply_style(); if they are missing, apply_style() raises an error
  rather than silently falling back to another font;
* PDF text kept as text: pdf.fonttype = ps.fonttype = 42 (TrueType embedded), svg.fonttype = 'none', so
  Adobe Illustrator opens labels as editable words; scatter clouds with many points are rasterised only where a
  script passes rasterized=True;
* colours: copied from the original notebooks
  - subpopulations / populations: sc_exploration/0_patient_stratification.ipynb (color_dict, col_lv2)
  - donor groups: same notebook (cluster_colors)
  - PD vs control: same notebook and the comparison plots of 3_match_modules_to_multiregion_pd.ipynb
    (PD #ffd652, control green); the Fullard UMAPs alone use #F2CB57 / #94D7F2 (CONDITION_COLORS_UMAP)
  - spatial (GloScope) donor groups: st_exploration/.../5_0_sample_rep_patpy.ipynb (cell 43)
  - Fullard regions: 3_match_modules_to_multiregion_pd.ipynb
* colormaps: diverging RdYlBu_r and sequential Purples, the two most used in the original notebooks.

Usage:
    from style import apply_style, save_figure, SUBTYPE_COLORS, GROUP_COLORS, CONDITION_COLORS, CMAP_DIV, CMAP_SEQ
    apply_style()                      # call after importing scanpy/seaborn (they can reset rcParams)
    save_figure(fig, FIG / "F01_name")  # writes F01_name.pdf (editable text) and F01_name.png
"""
from pathlib import Path
import matplotlib

FONT_FAMILY = "Arial"

from pdrep import config as cfg  # noqa: E402

FONT_DIR = cfg.LEGACY / "revision/env/fonts"  # Arial TTFs (not committed)
RC = {
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
    "font.family": "sans-serif", "font.sans-serif": [FONT_FAMILY],
    "mathtext.fontset": "custom", "mathtext.rm": FONT_FAMILY, "mathtext.it": FONT_FAMILY + ":italic",
    "mathtext.bf": FONT_FAMILY + ":bold",
    "axes.unicode_minus": True,
    "savefig.dpi": 300, "figure.dpi": 100,
    "savefig.facecolor": "white", "figure.facecolor": "white",
}

# ---- palettes (verbatim from the original notebooks) ------------------------------------------------------------
SUBTYPE_COLORS = {
    # astrocytes
    "A0": "#FADA7A", "AR1": "#F0A04B", "AR2": "#EB5A3C", "WMAs": "#6E5773", "ARPs": "#75a47f", "SVZCs": "#A4B465",
    # microglia / macrophages
    "MH": "#A8DADC", "MH_HSP": "#F4d58d", "DAM": "#B5A0c9", "M1": "#C7b198", "PM": "#DAA520",
    # oligodendroglia
    "MOL-A": "#B08BBB", "MOL-B": "#556FB5", "MOL-F": "#BEF0CB", "MOL-G": "#3A6351",
    "PDAO-1": "#804674", "PDAO-2": "#9ADCFF", "PDAO-3": "#76BA99",
    "immOL-A": "#E5E0FF", "immOL-B": "#B2A4FF", "immOL-C": "#767C77", "COPs": "#FFB4B4",
    "OPCs-A": "#D0B8A8", "OPCs-C": "#FFD966", "OPCs-D": "#CC7351", "PDAOPCs": "#FAAB78",
    # vascular, lymphocytes
    "Endothelial": "#4169E1", "Pericytes": "#C0C0C0", "VLMCs": "#20B2AA", "Fibroblasts": "#87CEFA",
    "B_cells": "#DC143C", "NK_cells": "#ADFF2F", "CD8_TRMM": "#6495ED", "CD8_TRME": "#1c55ba", "CD4_TM": "#BDB76B",
}
POPULATION_COLORS = {
    "Astrocytes": "#A8D8C3", "EpCs": "#5D866C", "Excitatory": "#E5073F", "Fibroblasts": "#C6DCF2",
    "Interneuron": "#FA5C3C", "Lymphocytes": "#AEE3E1", "MSN": "#b0f9a0", "Macrophages": "#fcaefa",
    "Microglia": "#9AA7E8", "OPCs": "#F2C49B", "Oligodendrocytes": "#F6E6A6", "Vascular": "#E7A1C4", "PM": "#D1BB9E",
}
GROUP_COLORS = {"Reactive": "#ae54a0", "UPR": "#429d78", "Other": "#fedbb6",
                "Mixed": "#8c8c8c", "Borderline": "#bdbdbd", "Unclassified": "#d9d9d9"}
CONDITION_COLORS = {"PD": "#ffd652", "Control": "#1ec167"}
CONDITION_COLORS_UMAP = {"PD": "#F2CB57", "Control": "#94D7F2"}     # Fullard UMAPs only
SPATIAL_GROUP_COLORS = {"Reactive": "#C04ABC", "UPR-high": "#419D78",
                        "Homeostatic (leiden 1)": "#E0A458", "Homeostatic (leiden 2)": "#FFDBB5"}
REGION_COLORS_FULLARD = {"medial globus pallidus": "#542b61", "primary motor cortex": "#973264",
                         "prefrontal cortex": "#cf5e56", "primary visual cortex": "#daab85",
                         "dorsal motor nucleus of vagus nerve": "#140a2c"}
COVARIATE_TYPE_COLORS = {"relevant": "#0173B2", "technical": "#DE8F05", "contextual": "#029E73"}  # seaborn 'colorblind', as in 5_0
BRAAK_COLORS = {"LBD Braak": "#4d4d4d", "NFT Braak": "#a6a6a6"}  # greys, as the Braak palette of the multiregion notebook
INK, MUTED, GRID, NULL_BAND = "#222222", "#6b6b6b", "#e6e6e6", "#c8c8c8"

# Outcomes and robustness (diverging around a neutral grey; validated: CVD dE >= 15, normal-vision dE >= 16.9).
OUTCOME_COLORS = {"reproduced": "#313695", "direction_only": "#5b8fc7", "not_reproduced": "#f4a259",
                  "contradicted": "#c0282d", "not_testable": "#efefef"}
ROBUSTNESS_COLORS = {"robust": "#313695", "supported": "#5b8fc7", "not_replicated": "#d9d9d9",
                     "inconsistent": "#f4a259", "contradicted": "#c0282d", "insufficient": "#efefef"}
OUTCOME_LABELS = {"reproduced": "reproduced", "direction_only": "direction only", "not_reproduced": "not reproduced",
                  "contradicted": "contradicted", "not_testable": "not testable"}
ROBUSTNESS_LABELS = {"robust": "robust", "supported": "supported", "not_replicated": "not replicated",
                     "inconsistent": "inconsistent", "contradicted": "contradicted", "insufficient": "insufficient data"}

# Claim types: dataviz reference categorical palette, slots 1-6 in fixed order (light-mode values).
CLAIM_GROUPS = ["composition", "gene expression", "cell state", "genetics", "communication / trajectory", "other"]
CLAIM_GROUP_COLORS = dict(zip(CLAIM_GROUPS, ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]))
TEMPLATE_GROUP = {"composition_donor": "composition", "de_pseudobulk": "gene expression", "geneset_shift": "gene expression",
                  "state_signature": "cell state", "genetic_enrichment": "genetics", "communication_liana":
                  "communication / trajectory", "trajectory": "communication / trajectory"}
REGION_GROUP = {"SN": "midbrain / SN", "midbrain": "midbrain / SN", "putamen": "striatum", "caudate": "striatum",
                "striatum": "striatum", "PFC": "cortex", "DLPFC": "cortex", "cingulate": "cortex", "temporal": "cortex",
                "MTG": "cortex", "motor": "cortex", "visual": "cortex", "amygdala": "other", "DMNX": "other",
                "GPi": "other", "cerebellum": "other"}

CMAP_DIV = "RdYlBu_r"     # diverging (z-scores, correlations)
CMAP_SEQ = "Purples"      # sequential (counts, fractions, R²)


def register_fonts():
    """Make the Arial TTFs in revision/env/fonts visible to matplotlib."""
    from matplotlib import font_manager
    files = sorted(FONT_DIR.glob("*.ttf")) + sorted(FONT_DIR.glob("*.TTF"))
    for f in files:
        font_manager.fontManager.addfont(str(f))
    names = {font_manager.FontProperties(fname=str(f)).get_name() for f in files}
    if FONT_FAMILY not in names:
        raise RuntimeError(f"{FONT_FAMILY} TTF files not found in {FONT_DIR}; copy Arial*.ttf there")
    return files


def apply_style(**overrides):
    """Set the manuscript rcParams. Extra rcParams (e.g. font sizes) can be passed as keyword arguments."""
    register_fonts()
    rc = dict(RC)
    rc.update({k.replace("__", "."): v for k, v in overrides.items()})
    matplotlib.rcParams.update(rc)
    from matplotlib import colorbar as _colorbar
    _colorbar.Colorbar.n_rasterize = 10 ** 9   # colorbars as vector cells, not bitmaps (Illustrator)
    return rc


def save_figure(fig, stem, png=True, **kw):
    """Save <stem>.pdf with editable text (and <stem>.png). Re-applies fonttype in case a library reset it."""
    matplotlib.rcParams.update({"pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none"})
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    kw.setdefault("bbox_inches", "tight")
    fig.savefig(stem.with_suffix(".pdf"), **kw)
    if png:
        fig.savefig(stem.with_suffix(".png"), dpi=kw.pop("dpi", 200), **kw)
    return stem.with_suffix(".pdf")


def heatmap(ax, M, cmap=None, vmin=None, vmax=None, origin="upper", aspect="auto"):
    """Heatmap drawn as vector cells (pcolormesh) instead of an embedded bitmap (imshow), so that it shows and stays
    editable in Adobe Illustrator. Cell (i, j) is centred on x = j, y = i, as with imshow."""
    import numpy as np
    M = np.ma.masked_invalid(np.ma.asarray(M, dtype=float))
    nr, nc = M.shape
    qm = ax.pcolormesh(np.arange(nc + 1) - 0.5, np.arange(nr + 1) - 0.5, M, cmap=cmap, vmin=vmin, vmax=vmax,
                       shading="flat", edgecolors="face", linewidth=0)
    ax.set_xlim(-0.5, nc - 0.5)
    ax.set_ylim((nr - 0.5, -0.5) if origin == "upper" else (-0.5, nr - 0.5))
    if aspect == "equal":
        ax.set_aspect("equal")
    return qm

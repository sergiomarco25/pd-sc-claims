"""The stress programme in tissue sections of Lewy-body dementia cortex (Visium; Jin et al. 2026).

Usage: python -I scripts/47_jin2026_visium.py
Data: GSE293896 (Jin et al. 2026, Acta Neuropathol), 10x Visium of temporal cortex, 6 Lewy-body dementia (two with SNCA
triplication) and 4 control brains, two batches. For each section the in-tissue spots are summed to one pseudobulk
profile and scored as in scripts/27_programmes.py (mean z-scored log2 CPM of the programme's genes across sections);
LBD sections are compared with controls (Hedges' g, Mann-Whitney). Per spot, the heat-shock score (mean log1p CP10k
of the ten shared heat-shock genes) is mapped. Writes results/biology/jin2026_visium.tsv and figures/B20_visium_cortex.
"""
from pathlib import Path
import importlib.util
import sys
import warnings
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.plot.style import apply_style, save_figure, CONDITION_COLORS, INK, MUTED  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

warnings.simplefilter("ignore")

from pdrep import config as cfg  # noqa: E402

RAW = cfg.RAW / "jin2026"
B = cfg.RESULTS / "biology"
spec = importlib.util.spec_from_file_location("prog", ROOT / "scripts" / "27_programmes.py")
prog = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prog)
HEAT10 = ["HSPA1A", "HSPA1B", "HSPH1", "DNAJB1", "DNAJA1", "HSP90AA1", "BAG3", "CHORDC1", "STIP1", "AHSA1"]
SHOW = ["shared heat-shock genes (10)", "heat response", "chaperone folding", "unfolded protein response", "hypoxia",
        "cold-inducible (CIRBP, RBM3)", "disease-associated oligodendrocytes (Adams)", "UPR astrocytes (Barba-Reyes)",
        "inflammatory response", "oxidative phosphorylation"]


def samples():
    rows = []
    for line in open(RAW / "samples.tsv"):
        f = line.rstrip("\n").split("\t")
        if len(f) < 4:
            continue
        rows.append({"gsm": f[0].split("= ")[1], "title": f[1].split("= ")[1], "batch": f[3].split(": ")[1]})
    s = pd.DataFrame(rows)
    s["lbd"] = s.title.str.startswith("LBD")
    s["label"] = s.title.str.replace(" brain", "").str.replace(" with ", " · ").str.replace(" and SNCA triplication", ", SNCA×3")
    return s


def load(gsm):
    pre = next(RAW.glob(f"{gsm}_*_matrix.mtx.gz")).name.replace("_matrix.mtx.gz", "")
    m = scipy.io.mmread(RAW / f"{pre}_matrix.mtx.gz").tocsc()  # genes x spots
    genes = pd.read_csv(RAW / f"{pre}_features.tsv.gz", sep="\t", header=None)[1].str.upper().values
    bc = pd.read_csv(RAW / f"{pre}_barcodes.tsv.gz", sep="\t", header=None)[0].values
    pos = pd.read_csv(RAW / f"{pre}_tissue_positions.csv.gz")
    if "barcode" not in pos.columns:
        pos.columns = ["barcode", "in_tissue", "row", "col", "y", "x"]
    pos = pos.set_index("barcode").reindex(bc)
    keep = pos.in_tissue.values == 1
    return m[:, keep], genes, pos[keep]


def main():
    s = samples()
    pbs, spots = {}, []
    for r in s.itertuples():
        m, genes, pos = load(r.gsm)
        pbs[r.gsm] = pd.Series(np.asarray(m.sum(1)).ravel(), index=genes).groupby(level=0).sum()
        cp = m.multiply(1e4 / np.maximum(np.asarray(m.sum(0)).ravel(), 1)).tocsr()
        idx = [i for i, g in enumerate(genes) if g in HEAT10]
        heat = np.asarray(np.log1p(cp[idx].toarray()).mean(0)).ravel()
        spots.append(pd.DataFrame({"gsm": r.gsm, "x": pos.iloc[:, -1].values, "y": pos.iloc[:, -2].values, "heat": heat}))
        print(r.gsm, m.shape, flush=True)
    pb = pd.DataFrame(pbs).fillna(0)
    lc = np.log2(pb.div(pb.sum(0), axis=1) * 1e6 + 1)
    lc = lc[(pb > 0).mean(1) == 1]
    z = lc.sub(lc.mean(1), axis=0).div(lc.std(1).replace(0, np.nan), axis=0)
    progs = {k: v[0] for k, v in prog.programmes().items()}
    progs["shared heat-shock genes (10)"] = HEAT10
    lbd = s.set_index("gsm").lbd.reindex(z.columns).values
    rows = []
    for p in SHOW:
        g = [k for k in progs[p] if k in z.index]
        sc_ = z.loc[g].mean(0)
        a, r_ = sc_[lbd], sc_[~lbd]
        spd = np.sqrt(((len(a) - 1) * a.var() + (len(r_) - 1) * r_.var()) / (len(a) + len(r_) - 2))
        rows.append({"programme": p, "n_lbd": len(a), "n_ctrl": len(r_), "hedges_g": (a.mean() - r_.mean()) / spd,
                     "p": stats.mannwhitneyu(a, r_).pvalue})
    t = pd.DataFrame(rows)
    t.to_csv(B / "jin2026_visium.tsv", sep="\t", index=False)
    print(t.round(3).to_string(index=False))
    figure(pd.concat(spots), s)


def figure(sp_, s):
    apply_style(font__size=7)
    order = list(s[~s.lbd].gsm) + list(s[s.lbd].gsm)
    hi = sp_.heat.quantile(0.99)
    fig, axs = plt.subplots(2, 5, figsize=(7.2, 3.4))
    for ax, g in zip(axs.ravel(), order):
        x = sp_[sp_.gsm == g]
        ax.scatter(x.x, -x.y, s=0.6, c=x.heat.clip(0, hi), cmap="magma_r", vmin=0, vmax=hi, lw=0, rasterized=True)
        ax.set_aspect("equal")
        ax.axis("off")
        r = s.set_index("gsm").loc[g]
        ax.set_title(f"{r.label}\nmedian {x.heat.median():.2f}", fontsize=5.5, loc="left",
                     color="#a50026" if r.lbd else "#118a50")
    sm = plt.cm.ScalarMappable(cmap="magma_r", norm=plt.Normalize(0, hi))
    cb = fig.colorbar(sm, ax=axs, fraction=0.015, pad=0.01)
    cb.set_label("heat-shock score per spot", fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    save_figure(fig, ROOT / "figures" / "B20_visium_cortex")
    plt.close(fig)


if __name__ == "__main__":
    main()

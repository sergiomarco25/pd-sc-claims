"""Is the stress programme higher in neurons that carry Lewy pathology? GeoMx whole-transcriptome data of Goralski et al.

Usage: python -I scripts/45_lewy_neurons_geomx.py
Data: Goralski et al. 2024 (Nat Commun), anterior cingulate cortex, NanoString GeoMx human whole-transcriptome atlas;
238 areas of illumination (AOIs) of neurons from layers 2/3, 5 and 6: NeuN-positive neurons without pSyn inclusions
in controls and cases (PD, PDD, DLB), and pSyn-positive neurons (Lewy pathology) in cases. Counts and annotations
were extracted from the deposited NanoStringGeoMxSet (raw/goralski2024/goralski_human_*.csv).
AOIs flagged in QC are dropped; probes are collapsed to genes; counts are normalised to library size (per 100,000),
log2-transformed and z-scored per gene across AOIs. Programme scores are the mean z of the programme's genes
(scripts/27_programmes.py gene sets, plus the ten shared heat-shock genes of this study), minus the mean score of 200
random gene sets matched gene by gene for mean expression (pSyn AOIs are much shallower than NeuN AOIs).
Tests: (1) Lewy-bearing vs Lewy-free neurons of the same case and layer (paired across case x layer units, Wilcoxon);
(2) Lewy-free neurons of cases vs control neurons (case means, Mann-Whitney), within layer and pooled.
Writes results/biology/lewy_neurons_geomx.tsv and figures/B19_lewy_neurons (PDF + PNG).
"""
from pathlib import Path
import importlib.util
import sys
import warnings
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.plot.style import apply_style, save_figure, CONDITION_COLORS, INK, MUTED, GRID  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

warnings.simplefilter("ignore")

from pdrep import config as cfg  # noqa: E402

RAW = cfg.RAW / "goralski2024"
B = cfg.RESULTS / "biology"
spec = importlib.util.spec_from_file_location("prog", ROOT / "scripts" / "27_programmes.py")
prog = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prog)
HEAT10 = ["HSPA1A", "HSPA1B", "HSPH1", "DNAJB1", "DNAJA1", "HSP90AA1", "BAG3", "CHORDC1", "STIP1", "AHSA1"]
SHOW = ["shared heat-shock genes (10)", "heat response", "chaperone folding", "unfolded protein response",
        "cold-inducible (CIRBP, RBM3)", "hypoxia", "oxidative phosphorylation", "inflammatory response"]


def load():
    e = pd.read_csv(RAW / "goralski_human_exprs.csv", index_col=0)
    ph = pd.read_csv(RAW / "goralski_human_pheno.csv", index_col=0)
    ft = pd.read_csv(RAW / "goralski_human_features.csv", index_col=0)
    keep = ph.QC.isna()
    e, ph = e.loc[:, keep], ph[keep]
    e = e[~ft.loc[e.index, "Negative"].astype(bool).values]
    e.index = ft.loc[e.index, "TargetName"].astype(str).str.upper().values
    e = e.groupby(level=0).sum()
    lx = np.log2(e / e.sum(0) * 1e5 + 1)  # library-size normalisation (Q3 is 0 in the shallow pSyn AOIs)
    lx = lx[(e >= 1).mean(1) >= 0.1]
    z = lx.sub(lx.mean(1), axis=0).div(lx.std(1).replace(0, np.nan), axis=0).dropna()
    ph["group"] = np.where(ph.segment.eq("pSyn"), "Lewy-bearing neurons",
                           np.where(ph.Diagnosis.eq("Ctrl"), "control neurons", "Lewy-free neurons (cases)"))
    return z, ph, lx.mean(1).loc[z.index]


def main():
    z, ph, gmean = load()
    progs = {k: v[0] for k, v in prog.programmes().items()}
    progs["shared heat-shock genes (10)"] = HEAT10
    # pSyn AOIs hold ~5 nuclei and ~10x fewer counts than NeuN AOIs, which shifts z-scores of all genes; each
    # programme score is therefore expressed relative to 200 random gene sets of the same size, matched gene by gene
    # for mean expression (deciles), scored in the same AOIs
    rng = np.random.default_rng(0)
    dec = pd.qcut(gmean, 10, labels=False)
    pools = {k: list(v.index) for k, v in dec.groupby(dec)}
    sc = {}
    for k, genes in progs.items():
        g = [x for x in genes if x in z.index]
        if len(g) < (2 if len(genes) == 2 else 5):
            continue
        obs_ = z.loc[g].mean(0)
        rnd = np.mean([z.loc[[rng.choice(pools[dec[x]]) for x in g]].mean(0).values for _ in range(200)], axis=0)
        sc[k] = obs_ - rnd
    sc = pd.DataFrame(sc)
    d = ph.join(sc)
    rows = []
    for p in SHOW:
        # (1) paired within case x layer
        u = d.groupby(["Case", "Layer", "group"])[p].mean().unstack()
        pr = u.dropna(subset=["Lewy-bearing neurons", "Lewy-free neurons (cases)"])
        diff = pr["Lewy-bearing neurons"] - pr["Lewy-free neurons (cases)"]
        w = stats.wilcoxon(diff) if len(diff) >= 5 else None
        rows.append({"programme": p, "test": "Lewy-bearing vs Lewy-free neurons (same case and layer)",
                     "n": len(diff), "n_cases": pr.index.get_level_values(0).nunique(), "mean_diff_sd": diff.mean(),
                     "p": w.pvalue if w else np.nan})
        # (2) Lewy-free neurons of cases vs control neurons, case means
        cm = d[d.group != "Lewy-bearing neurons"].groupby(["Case", "group"])[p].mean().reset_index()
        a = cm[cm.group.eq("Lewy-free neurons (cases)")][p]
        r = cm[cm.group.eq("control neurons")][p]
        mw = stats.mannwhitneyu(a, r)
        rows.append({"programme": p, "test": "Lewy-free neurons of cases vs control neurons", "n": len(a) + len(r),
                     "n_cases": len(a), "mean_diff_sd": a.mean() - r.mean(), "p": mw.pvalue})
    t = pd.DataFrame(rows)
    t.to_csv(B / "lewy_neurons_geomx.tsv", sep="\t", index=False)
    d.to_csv(B / "lewy_neurons_geomx_aoi.tsv", sep="\t")
    pd.set_option("display.width", 250)
    print(d.groupby(["Diagnosis", "segment"]).size().to_string())
    print(t.round(4).to_string(index=False))
    figure(d, t)


def figure(d, t):
    apply_style(font__size=7, axes__linewidth=0.6, axes__edgecolor=INK)
    show = ["shared heat-shock genes (10)", "chaperone folding", "unfolded protein response", "cold-inducible (CIRBP, RBM3)"]
    cols = {"control neurons": CONDITION_COLORS["Control"], "Lewy-free neurons (cases)": "#d4a400",
            "Lewy-bearing neurons": "#a50026"}
    fig, axs = plt.subplots(1, len(show), figsize=(7.2, 2.9), gridspec_kw={"wspace": 0.45})
    order = list(cols)
    for ax, p in zip(axs, show):
        u = d.groupby(["Case", "Layer", "group"])[p].mean().reset_index()
        for i, g in enumerate(order):
            v = u[u.group == g][p]
            ax.scatter(np.full(len(v), i) + np.random.default_rng(i).uniform(-0.15, 0.15, len(v)), v, s=8, color=cols[g],
                       edgecolor=INK, lw=0.2, zorder=3)
            ax.plot([i - 0.25, i + 0.25], [v.median()] * 2, color=INK, lw=1.2, zorder=4)
        pr = u.pivot_table(index=["Case", "Layer"], columns="group", values=p)
        pr = pr.dropna(subset=[c for c in ["Lewy-bearing neurons", "Lewy-free neurons (cases)"] if c in pr])
        for _, r in pr.iterrows():  # lines join the two neuron types of the same case and layer
            ax.plot([1, 2], [r["Lewy-free neurons (cases)"], r["Lewy-bearing neurons"]], color=GRID, lw=0.6, zorder=1)
        tt = t.set_index(["programme", "test"])
        p1 = tt.loc[(p, "Lewy-bearing vs Lewy-free neurons (same case and layer)"), "p"]
        p2 = tt.loc[(p, "Lewy-free neurons of cases vs control neurons"), "p"]
        ax.set_title(f"{p}\nLewy+ vs Lewy−: p = {p1:.2g}\ncases vs controls: p = {p2:.2g}", fontsize=6, loc="left")
        ax.set_xticks(range(3), ["control", "case,\nLewy−", "case,\nLewy+"], fontsize=6)
        for s_ in ["top", "right"]:
            ax.spines[s_].set_visible(False)
        ax.tick_params(labelsize=5.5, length=2)
    axs[0].set_ylabel("programme score (z)\nper case × layer", fontsize=6)
    fig.subplots_adjust(top=0.75, bottom=0.15, left=0.08, right=0.99)
    save_figure(fig, ROOT / "figures" / "B19_lewy_neurons")
    plt.close(fig)


if __name__ == "__main__":
    main()

"""Is the PD stress state specific to Lewy-body disease? Programme scores across all PsychAD diagnoses (DLPFC).

Usage: python -I scripts/41b_psychad_specificity.py
Input: results/psychad/pb__<class>.csv.gz, donors.tsv and qc.tsv (scripts/41a_psychad_pseudobulk.py). Donor units with
at least 20 nuclei of a class are scored as in scripts/27_programmes.py (mean z-scored log2 CPM of the programme's genes
across all donors of the class). For every diagnostic group versus neurotypical controls, the programme score is
modelled as score ~ group + age + sex + log median genes per nucleus (ordinary least squares); the group coefficient
is expressed in control standard deviations, with its standard error and p-value. Also reported: the share of the
group's donors above the controls' 95th percentile (raw scores). A donor stress index averages the heat-response and
chaperone-folding scores of astrocytes, microglia, oligodendrocytes and OPCs (each in control SD units).
Writes results/psychad/specificity.tsv and figures/B17_specificity (PDF + PNG).
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
from pdrep.plot.style import apply_style, save_figure, heatmap, CMAP_DIV, INK, MUTED, GRID  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

warnings.simplefilter("ignore")

from pdrep import config as cfg  # noqa: E402

P = cfg.RESULTS / "psychad"
spec = importlib.util.spec_from_file_location("prog", ROOT / "scripts" / "27_programmes.py")
prog = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prog)
GROUPS = ["Lewy-body disease", "Lewy-body disease with AD", "Alzheimer's disease", "vascular dementia", "tauopathy",
          "frontotemporal dementia", "other dementia", "schizophrenia with dementia", "schizophrenia"]
ROWS = [("heat response", "Oligodendrocytes"), ("heat response", "Astrocytes"), ("heat response", "Microglia"),
        ("chaperone folding", "Oligodendrocytes"), ("chaperone folding", "Microglia"),
        ("unfolded protein response", "Oligodendrocytes"), ("hypoxia", "Microglia"),
        ("cold-inducible (CIRBP, RBM3)", "Oligodendrocytes"), ("cold-inducible (CIRBP, RBM3)", "Astrocytes"),
        ("UPR astrocytes (Barba-Reyes)", "Astrocytes"), ("UPR astrocytes, Astrocytes3 (Martirosyan)", "Astrocytes"),
        ("UPR microglia (Barba-Reyes)", "Microglia"), ("disease-associated oligodendrocytes (Adams)", "Oligodendrocytes"),
        ("GPNMB-high microglia (Smajić)", "Microglia"), ("AR2 astrocytes (Barba-Reyes)", "Astrocytes"),
        ("inflammatory response", "Microglia"), ("heat response", "Excitatory neurons"),
        ("oxidative phosphorylation", "Excitatory neurons"), ("stress index", "glia")]


def scores():
    progs = prog.programmes()
    qc = pd.read_csv(P / "qc.tsv", sep="\t")
    out = []
    for f in sorted(P.glob("pb__*.csv.gz")):
        cls = f.name[len("pb__"):-len(".csv.gz")].replace("_", " ")
        x = pd.read_csv(f, index_col=0)
        x.index = x.index.str.upper()
        q = qc[qc.cell_class.eq(cls)].set_index("donor")
        keep = [d for d in x.columns if q.nuclei.get(d, 0) >= 20]
        x = x[keep]
        lc = np.log2(x.div(x.sum(0), axis=1) * 1e6 + 1)
        lc = lc[(x > 0).mean(1) >= 0.1]
        z = lc.sub(lc.mean(1), axis=0).div(lc.std(1).replace(0, np.nan), axis=0)
        for pname, (genes, classes) in progs.items():
            if cls not in classes:
                continue
            g = [k for k in genes if k in z.index]
            if len(g) < (2 if len(genes) == 2 else 5):
                continue
            out.append(pd.DataFrame({"donor": z.columns, "cell_class": cls, "programme": pname,
                                     "score": z.loc[g].mean(0).values,
                                     "lg": np.log(q.loc[z.columns, "median_genes"].values)}))
    return pd.concat(out, ignore_index=True)


def fit(d, grp):
    x = d[d.group.isin([grp, "control"])].dropna(subset=["score", "age", "lg"])
    if (x.group == grp).sum() < 5:
        return None
    ctrl = x[x.group == "control"].score
    X = np.c_[np.ones(len(x)), (x.group == grp).astype(float), x.age, (x.sex == "male").astype(float), x.lg]
    y = x.score.values
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ beta
    dof = len(y) - X.shape[1]
    s2 = res @ res / dof
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    t = beta[1] / se[1]
    sd = ctrl.std(ddof=1)
    above = (x[x.group == grp].score > np.quantile(ctrl, 0.95)).mean()
    return {"effect_sd": beta[1] / sd, "se_sd": se[1] / sd, "p": 2 * stats.t.sf(abs(t), dof),
            "n_group": int((x.group == grp).sum()), "n_ctrl": len(ctrl), "share_above_ctrl95": above}


def main():
    meta = pd.read_csv(P / "donors.tsv", sep="\t")
    s = scores().merge(meta[["donor", "group", "age", "sex"]], on="donor")
    # donor stress index: heat response + chaperone folding of four glial classes, in control SD units
    gl = s[s.programme.isin(["heat response", "chaperone folding"])
           & s.cell_class.isin(["Astrocytes", "Microglia", "Oligodendrocytes", "OPCs"])].copy()
    ref = gl[gl.group == "control"].groupby(["programme", "cell_class"]).score.agg(["mean", "std"])
    gl = gl.join(ref, on=["programme", "cell_class"])
    gl["zs"] = (gl.score - gl["mean"]) / gl["std"]
    idx = gl.groupby("donor").agg(score=("zs", "mean"), lg=("lg", "mean")).reset_index()
    idx = idx.merge(meta[["donor", "group", "age", "sex"]], on="donor").assign(programme="stress index", cell_class="glia")
    s = pd.concat([s, idx], ignore_index=True)
    rows = []
    for (p, c) in ROWS:
        d = s[(s.programme == p) & (s.cell_class == c)]
        for g in GROUPS:
            r = fit(d, g)
            if r:
                rows.append({"programme": p, "cell_class": c, "group": g, **r})
    t = pd.DataFrame(rows)
    t["fdr"] = stats.false_discovery_control(t.p)
    t.to_csv(P / "specificity.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(t.pivot_table(index=["programme", "cell_class"], columns="group", values="effect_sd").round(2)
          .reindex(columns=GROUPS).to_string())
    figure(t, s)


def figure(t, s):
    apply_style(font__size=7, axes__linewidth=0.6, axes__edgecolor=INK)
    n = t.drop_duplicates("group").set_index("group").n_group
    fig = plt.figure(figsize=(7.2, 5.6))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.6, 1], wspace=0.75)
    ax = fig.add_subplot(gs[0])
    M = t.pivot_table(index=["programme", "cell_class"], columns="group", values="effect_sd").reindex(columns=GROUPS)
    F = t.pivot_table(index=["programme", "cell_class"], columns="group", values="fdr").reindex(columns=GROUPS)
    M, F = M.reindex(ROWS), F.reindex(ROWS)
    qm = heatmap(ax, M.values, cmap=CMAP_DIV, vmin=-1, vmax=1)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if F.values[i, j] < 0.05:
                ax.text(j, i, "*", ha="center", va="center", fontsize=8, color="white" if abs(M.values[i, j]) > 0.6 else INK)
    ax.set_xticks(range(len(GROUPS)), [f"{g} ({n.get(g, 0)})" for g in GROUPS], rotation=60, ha="right", fontsize=6)
    lab = [f"{p} · {c.lower()}" if c != "glia" else "stress index · glia" for p, c in ROWS]
    ax.set_yticks(range(len(ROWS)), lab, fontsize=6)
    ax.axvline(1.5, color="white", lw=2)
    ax.axhline(len(ROWS) - 1.5, color="white", lw=2)
    cb = fig.colorbar(qm, ax=ax, fraction=0.025, pad=0.33, orientation="horizontal", aspect=40)
    cb.set_label("difference from controls (control SD), adjusted for age, sex and nucleus quality", fontsize=6)
    cb.ax.tick_params(labelsize=6)
    ax.set_title("a  Programmes across diagnoses (PsychAD, DLPFC)", loc="left", fontsize=8, fontweight="bold")
    ax = fig.add_subplot(gs[1])
    d = s[(s.programme == "stress index")]
    order = ["control"] + GROUPS
    for i, g in enumerate(order):
        v = d[d.group == g].score.dropna()
        if not len(v):
            continue
        jit = np.random.default_rng(i).uniform(-0.25, 0.25, len(v))
        col = "#1ec167" if g == "control" else ("#d4a400" if g.startswith("Lewy") else "#8c8c8c")
        ax.scatter(v.clip(-3, 6), i + jit, s=3, color=col, alpha=0.7, lw=0, rasterized=True)
        ax.plot([v.median()] * 2, [i - 0.35, i + 0.35], color=INK, lw=1.2)
    ax.axvline(0, color=GRID, lw=0.8)
    ax.set_yticks(range(len(order)), order, fontsize=6)
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_xlabel("donor stress index (control SD)")
    for s_ in ["top", "right"]:
        ax.spines[s_].set_visible(False)
    ax.set_title("b  Donor stress index", loc="left", fontsize=8, fontweight="bold")
    ax.set_xlim(-3, 6)
    fig.subplots_adjust(left=0.25, right=0.98, bottom=0.12, top=0.94)
    save_figure(fig, ROOT / "figures" / "B17_specificity")
    plt.close(fig)


if __name__ == "__main__":
    main()

"""Side analysis: are the two replicated patient axes also axes of normal variation, and do they track clinical course?

Usage: python -I scripts/subgroups/03_axes_controls_clinical.py [--figs]  (--figs: redraw figures from saved tables)
Builds the donor profiles of scripts/subgroups/01_patient_subgroups.py for ALL donors (cases and controls; every
feature in standard deviations from the controls of the same dataset x region, so controls are centred at 0).
Every donor of every dataset is scored on the two consensus axes: the mean PC1 and PC2 loadings of the four large
cohorts (results/subgroups/axis_loadings.tsv), applied to the raw control-referenced features (missing features
skipped, weights renormalised; donors with < 60% of the loading weight available are dropped).
  1. Normal variation: PCA of controls alone in each large cohort; congruence of control PC1-3 with the case axes.
  2. Shift: cases vs controls on each axis per dataset (Hedges' g), pooled across datasets (DerSimonian-Laird).
  3. Clinical course, within cases: disease duration and age at onset (partial Spearman given age at death), Lewy
     stage (Spearman), dementia (PDD/DLB) vs PD (rank-biserial r), per dataset with >= 10 informative cases,
     pooled across datasets (Fisher z, DerSimonian-Laird).
Writes results/subgroups/axis_scores_all.tsv, axes_controls.tsv, axes_shift.tsv, axes_clinical.tsv and
figures/subgroups/S3_axes_patients_controls, S4_axis_definitions, S5_axes_clinical.
"""
from pathlib import Path
import importlib.util
import sys
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.meta import dersimonian_laird  # noqa: E402
from pdrep.plot.style import apply_style, save_figure, CONDITION_COLORS, INK, MUTED, GRID  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

warnings.simplefilter("ignore")
spec = importlib.util.spec_from_file_location("sg", ROOT / "scripts" / "subgroups" / "01_patient_subgroups.py")
sg = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sg)

from pdrep import config as cfg  # noqa: E402

R = cfg.RESULTS
OUT, FIG = R / "subgroups", ROOT / "figures" / "subgroups"
META = cfg.LEGACY / "Human_Samples_Combined_Metadata_v3.xlsx"
LARGE = ["lee2026", "menon2026", "prashant2024", "barbareyes2026"]
CASES = {"PD", "PDD", "DLB", "LBD"}
CASE_C, CTRL_C = "#d4a400", CONDITION_COLORS["Control"]
AXIS_NAME = {"PC1": "axis 1: overall stress level", "PC2": "axis 2: proteostatic stress ↔ inflammation"}
LABEL = {"lee2026": "Lee 2026 (DLPFC)", "menon2026": "Menon 2026 (5 regions)", "prashant2024": "Prashant 2024 (5 regions)",
         "barbareyes2026": "Barba-Reyes 2026 (striatum)", "feleke2021": "Feleke 2021", "wang2024": "Wang 2024",
         "volpato2025": "Volpato 2025", "smajic2022": "Smajić 2022", "ma2025": "Ma 2025", "chia2025": "Chia 2025"}


def all_features():
    sg.CASES = CASES | {"control"}  # keep controls in the donor profiles
    X, meta = sg.features()
    return X, meta


def axis_scores(X, w):
    rows = []
    for ds in X.index.get_level_values(0).unique():
        m = X.loc[ds]
        for ax in ["PC1", "PC2"]:
            ww = w[ax].reindex(m.columns).dropna()
            M = m[ww.index]
            avail = M.notna().mul(ww ** 2, axis=1).sum(1) / (w[ax] ** 2).sum()
            sc = M.fillna(0).mul(ww, axis=1).sum(1) / np.sqrt(M.notna().mul(ww ** 2, axis=1).sum(1))
            rows.append(pd.DataFrame({"dataset": ds, "donor": m.index, "axis": ax, "score": sc.values,
                                      "weight_available": avail.values}))
    s = pd.concat(rows)
    return s[s.weight_available >= 0.6]


def shift_table(s):
    """Axis scores with a case/control group -> cases vs controls per dataset and axis (Hedges' g; >= 5 cases and
    >= 4 controls), pooled per axis across datasets (DerSimonian-Laird)."""
    sh = []
    for (ds, ax), d in s.groupby(["dataset", "axis"]):
        a, r = d[d.group == "case"].score, d[d.group == "control"].score
        if len(a) >= 5 and len(r) >= 4:
            sp = np.sqrt(((len(a) - 1) * a.var() + (len(r) - 1) * r.var()) / (len(a) + len(r) - 2))
            g = (a.mean() - r.mean()) / sp
            sh.append({"dataset": ds, "axis": ax, "g": g, "se": np.sqrt(1 / len(a) + 1 / len(r) + g ** 2 / (2 * (len(a) + len(r)))),
                       "n_case": len(a), "n_ctrl": len(r), "sd_ratio": a.std() / r.std()})
    sh = pd.DataFrame(sh)
    pooled = [{"dataset": "pooled", "axis": ax, **{k: v for k, v in dersimonian_laird(d.g, d.se).items()}}
              for ax, d in sh.groupby("axis")]
    return pd.concat([sh, pd.DataFrame(pooled)], ignore_index=True)


def clinical(meta):
    dur = pd.read_csv(R / "survival" / "durations_used.tsv", sep="\t", dtype={"donor": str})
    ons = []
    x = pd.read_excel(META)
    ons.append(pd.DataFrame({"dataset": "barbareyes2026", "donor": x["Subject ID"].astype(str),
                             "onset": pd.to_numeric(x["Age of onset"], errors="coerce")}))
    e = pd.read_csv(R / "survival" / "duration_extra.tsv", sep="\t", dtype={"donor": str})
    ons.append(e[["dataset", "donor"]].assign(onset=pd.to_numeric(e.age_onset, errors="coerce")))
    ons = pd.concat(ons).dropna()
    c = meta.reset_index().merge(dur, on=["dataset", "donor"], how="left").merge(ons, on=["dataset", "donor"], how="left")
    c.loc[c.dataset.eq("menon2026") & c.donor.eq("BANN1644"), "onset"] = np.nan
    return c


def partial_r(y, x, z):
    ok = np.isfinite(y) & np.isfinite(x) & np.isfinite(z)
    if ok.sum() < 10:
        return np.nan, int(ok.sum())
    ry, rx, rz = (stats.rankdata(v[ok]) for v in (y, x, z))
    res = lambda v: v - np.polyval(np.polyfit(rz, v, 1), rz)
    return stats.pearsonr(res(ry), res(rx))[0], int(ok.sum())


def pool_r(df):
    out = []
    for (ax, var), d in df.dropna(subset=["r"]).groupby(["axis", "variable"]):
        z = np.arctanh(d.r.clip(-0.99, 0.99))
        se = 1 / np.sqrt(d.n - 3)
        r = dersimonian_laird(z, se) if len(d) > 1 else {"est": z.iloc[0], "se": se.iloc[0],
                                                            "p": 2 * stats.norm.sf(abs(z.iloc[0] / se.iloc[0])), "I2": np.nan}
        out.append({"axis": ax, "variable": var, "dataset": "pooled", "r": np.tanh(r["est"]),
                    "lo": np.tanh(r["est"] - 1.96 * r["se"]), "hi": np.tanh(r["est"] + 1.96 * r["se"]),
                    "p": r["p"], "n": int(d.n.sum()), "k": len(d), "I2": r.get("I2", np.nan)})
    return pd.DataFrame(out)


def main():
    apply_style(font__size=7, axes__linewidth=0.6, axes__edgecolor=INK)
    X, meta = all_features()
    load = pd.read_csv(OUT / "axis_loadings.tsv", sep="\t", header=[0, 1], index_col=0)
    w = pd.DataFrame({"PC1": load[("mean", "PC1")], "PC2": load[("mean", "PC2")]})
    s = axis_scores(X, w).merge(meta.reset_index()[["dataset", "donor", "disease", "sex", "age", "stage"]],
                                on=["dataset", "donor"])
    s["group"] = np.where(s.disease.isin(CASES), "case", np.where(s.disease.eq("control"), "control", "other"))
    s.to_csv(OUT / "axis_scores_all.tsv", sep="\t", index=False)
    # 1. axes in controls alone
    ctl = []
    for ds in LARGE:
        m = X.loc[ds]
        m = m[meta.loc[ds].reindex(m.index).disease.eq("control").values]
        shared = [f for f in w.index if f in m.columns and m[f].notna().mean() > 0.7]
        Z = m[shared].fillna(0)
        Z = (Z - Z.mean()) / Z.std(ddof=1).replace(0, 1)
        p = PCA(3).fit(Z)
        for i in range(3):
            for ax in ["PC1", "PC2"]:
                ww = w.loc[shared, ax].values
                ctl.append({"cohort": ds, "controls": len(Z), "control_pc": i + 1, "var_explained": p.explained_variance_ratio_[i],
                            "case_axis": ax, "congruence": abs(sg.tucker(p.components_[i], ww))})
    ctl = pd.DataFrame(ctl)
    ctl.to_csv(OUT / "axes_controls.tsv", sep="\t", index=False)
    print(ctl.pivot_table(index=["cohort", "control_pc"], columns="case_axis", values="congruence").round(2).to_string())
    # 2. shift cases vs controls
    sh = shift_table(s)
    sh.to_csv(OUT / "axes_shift.tsv", sep="\t", index=False)
    print(sh.round(3).to_string(index=False))
    # 3. clinical, within cases
    c = clinical(meta)
    sc = s[s.group == "case"].merge(c[["dataset", "donor", "duration", "onset"]], on=["dataset", "donor"], how="left")
    rows = []
    for (ds, ax), d in sc.groupby(["dataset", "axis"]):
        age = pd.to_numeric(d.age, errors="coerce").values
        y = d.score.values
        for var, v in [("disease duration", d.duration.values), ("age at onset", d.onset.values)]:
            r, n = partial_r(y, v.astype(float), age)
            if np.isfinite(r):
                rows.append({"dataset": ds, "axis": ax, "variable": var, "r": r, "n": n})
        st = pd.to_numeric(d.stage, errors="coerce").values
        ok = np.isfinite(st)
        if ok.sum() >= 10 and len(np.unique(st[ok])) > 1:
            rows.append({"dataset": ds, "axis": ax, "variable": "Lewy stage", "r": stats.spearmanr(y[ok], st[ok])[0],
                         "n": int(ok.sum())})
        dem = d.disease.isin(["PDD", "DLB"]).values
        if dem.sum() >= 5 and (~dem).sum() >= 5:
            u = stats.mannwhitneyu(y[dem], y[~dem]).statistic
            rows.append({"dataset": ds, "axis": ax, "variable": "dementia vs PD",
                         "r": 2 * u / (dem.sum() * (~dem).sum()) - 1, "n": len(y)})
    cl = pd.DataFrame(rows)
    cl["se"] = 1 / np.sqrt(cl.n - 3)
    cl["lo"], cl["hi"] = np.tanh(np.arctanh(cl.r) - 1.96 * cl.se), np.tanh(np.arctanh(cl.r) + 1.96 * cl.se)
    cl["p"] = 2 * stats.norm.sf(abs(np.arctanh(cl.r) / cl.se))
    allc = pd.concat([cl, pool_r(cl)], ignore_index=True)
    allc.to_csv(OUT / "axes_clinical.tsv", sep="\t", index=False)
    print(allc[["axis", "variable", "dataset", "r", "p", "n"]].round(3).to_string(index=False))
    figures()


def figures():
    apply_style(font__size=7, axes__linewidth=0.6, axes__edgecolor=INK)
    load = pd.read_csv(OUT / "axis_loadings.tsv", sep="\t", header=[0, 1], index_col=0)
    s = pd.read_csv(OUT / "axis_scores_all.tsv", sep="\t", dtype={"donor": str})
    fig_axes(s, pd.read_csv(OUT / "axes_shift.tsv", sep="\t"))
    fig_definitions(load, pd.read_csv(OUT / "axes_controls.tsv", sep="\t"))
    fig_clinical(pd.read_csv(OUT / "axes_clinical.tsv", sep="\t"))


def fig_axes(s, sh):
    fig = plt.figure(figsize=(7.2, 5.0))
    gs = fig.add_gridspec(2, 4, height_ratios=[1.15, 1], hspace=0.55, wspace=0.35)
    for j, ds in enumerate(LARGE):
        ax = fig.add_subplot(gs[0, j])
        d = s[s.dataset == ds].pivot_table(index=["donor", "group"], columns="axis", values="score").reset_index()
        for g, col, z in [("control", CTRL_C, 2), ("case", CASE_C, 3)]:
            x = d[d.group == g]
            ax.scatter(x.PC1, x.PC2, s=7, color=col, edgecolor=INK, lw=0.2, alpha=0.85, zorder=z,
                       label="controls" if g == "control" else "Lewy-body cases")
        ax.axhline(0, color=GRID, lw=0.6, zorder=0)
        ax.axvline(0, color=GRID, lw=0.6, zorder=0)
        nc, nk = (d.group == "case").sum(), (d.group == "control").sum()
        ax.set_title(f"{LABEL[ds]}\n{nc} cases, {nk} controls", fontsize=6.3, loc="left")
        ax.set_xlabel("axis 1: overall stress", fontsize=6)
        if j == 0:
            ax.set_ylabel("axis 2: inflammation ← → proteostatic stress", fontsize=6)
        ax.tick_params(labelsize=5.5, length=2)
        if j == 0:
            fig.legend(*ax.get_legend_handles_labels(), frameon=False, fontsize=6.5, loc="upper right", ncol=2,
                       bbox_to_anchor=(0.99, 1.0), markerscale=1.8)
        for sp_ in ["top", "right"]:
            ax.spines[sp_].set_visible(False)
    for i, axn in enumerate(["PC1", "PC2"]):
        ax = fig.add_subplot(gs[1, 2 * i:2 * i + 2])
        d = sh[(sh.axis == axn)].copy()
        d = d[d.dataset != "pooled"].sort_values("n_case")
        y = np.arange(len(d))
        ax.errorbar(d.g, y, xerr=1.96 * d.se, fmt="o", color=INK, ms=3, lw=0.8)
        pdd = sh[(sh.axis == axn) & (sh.dataset == "pooled")].iloc[0]
        ax.errorbar(pdd.est, len(d) + 0.6, xerr=1.96 * pdd.se, fmt="D", color="#a50026", ms=4.5, lw=1.4)
        ax.set_yticks(list(y) + [len(d) + 0.6], [f"{x} ({int(a)}/{int(b)})" for x, a, b in zip(d.dataset, d.n_case, d.n_ctrl)]
                      + [f"pooled (p = {pdd.p:.1g})"], fontsize=5.5)
        ax.axvline(0, color=GRID, lw=0.8, zorder=0)
        ax.set_xlabel("cases vs controls (Hedges' g)", fontsize=6)
        ax.set_title(AXIS_NAME[axn], fontsize=6.5, loc="left")
        ax.tick_params(labelsize=5.5, length=2)
        for sp_ in ["top", "right"]:
            ax.spines[sp_].set_visible(False)
    fig.text(0.01, 0.985, "a  Patients and controls on the two axes (0 = average control)", fontsize=8, fontweight="bold", va="top")
    fig.text(0.01, 0.47, "b  How far cases sit from controls on each axis, per dataset (cases/controls)", fontsize=8,
             fontweight="bold", va="top")
    fig.subplots_adjust(top=0.875, bottom=0.09, left=0.17, right=0.98)
    save_figure(fig, FIG / "S3_axes_patients_controls")
    plt.close(fig)


def fig_definitions(load, ctl):
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[3.2, 1.1], hspace=0.32, wspace=1.45)
    axs = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]
    cohorts = [c for c in LARGE]
    for ax, axn in zip(axs, ["PC1", "PC2"]):
        L = pd.DataFrame({c: load[(c, axn)] for c in cohorts})
        m = load[("mean", axn)]
        top = pd.concat([m.sort_values().head(9), m.sort_values().tail(9)])
        y = np.arange(len(top))
        cols = ["#4575b4" if v < 0 else "#a50026" for v in top]
        ax.barh(y, top.values, color=cols, height=0.65)
        for k, f in enumerate(top.index):  # per-cohort loadings as dots
            ax.scatter(L.loc[f].values, np.full(len(cohorts), k), s=5, color=INK, zorder=3, lw=0)
        lab = [f.replace("disease-associated oligodendrocytes", "DA oligodendrocytes").replace("share · ", "share of ")
               .replace("_total", "") for f in top.index]
        ax.set_yticks(y, lab, fontsize=5.6)
        ax.axvline(0, color=INK, lw=0.6)
        ax.set_xlabel("loading (bar: mean; dots: each cohort)", fontsize=6)
        ax.set_title(("a  " if axn == "PC1" else "b  ") + AXIS_NAME[axn], fontsize=7, loc="left", fontweight="bold")
        ax.tick_params(labelsize=5.5, length=2)
        for sp_ in ["top", "right"]:
            ax.spines[sp_].set_visible(False)
    ax = fig.add_subplot(gs[1, :])
    best = ctl.sort_values("congruence").groupby(["cohort", "case_axis"]).tail(1)
    xs = np.arange(len(LARGE))
    for k, (axn, col) in enumerate([("PC1", "#a50026"), ("PC2", "#4575b4")]):
        b = best[best.case_axis == axn].set_index("cohort").reindex(LARGE)
        ax.bar(xs + (k - 0.5) * 0.36, b.congruence, width=0.34, color=col,
               label=AXIS_NAME[axn].split(":")[0] + " found in controls")
        for x0, r in zip(xs, b.itertuples()):
            ax.text(x0 + (k - 0.5) * 0.36, r.congruence + 0.02, f"control PC{r.control_pc}\n{r.var_explained:.0%} var.",
                    ha="center", va="bottom", fontsize=4.8)
    ax.axhline(0.85, color=MUTED, lw=0.6, ls="--")
    ax.text(len(LARGE) - 0.45, 0.86, "0.85 = same axis", fontsize=5, color=MUTED, ha="right", va="bottom")
    ax.set_xticks(xs, [f"{LABEL[c]}\n{int(ctl[ctl.cohort == c].controls.iloc[0])} controls" for c in LARGE], fontsize=5.6)
    ax.set_ylim(0, 1.3)
    ax.set_yticks([0, 0.5, 1])
    ax.set_ylabel("congruence with the case axis", fontsize=6)
    ax.set_title("c  PCA of controls alone: do the same axes appear without disease?", fontsize=7, loc="left",
                 fontweight="bold")
    ax.legend(frameon=False, fontsize=5.5, loc="lower right", ncol=2, bbox_to_anchor=(1.0, 1.08))
    ax.tick_params(labelsize=5.5, length=2)
    for sp_ in ["top", "right"]:
        ax.spines[sp_].set_visible(False)
    fig.subplots_adjust(left=0.24, right=0.98, top=0.95, bottom=0.08)
    save_figure(fig, FIG / "S4_axis_definitions")
    plt.close(fig)


def fig_clinical(allc):
    vars_ = ["disease duration", "age at onset", "Lewy stage", "dementia vs PD"]
    fig, axs = plt.subplots(1, 2, figsize=(7.2, 3.8), gridspec_kw={"wspace": 0.9})
    for ax, axn in zip(axs, ["PC2", "PC1"]):
        d = allc[allc.axis == axn]
        ticks, labs, y = [], [], 0
        for v in vars_:
            dv = d[d.variable == v]
            if dv.empty:
                continue
            for r in dv[dv.dataset != "pooled"].itertuples():
                ax.errorbar(r.r, y, xerr=[[r.r - r.lo], [r.hi - r.r]], fmt="o", color=MUTED, ms=3, lw=0.7)
                ticks.append(y)
                labs.append(f"   {r.dataset} (n = {int(r.n)})")
                y += 1
            pr = dv[dv.dataset == "pooled"]
            if len(pr):
                pr = pr.iloc[0]
                col = "#a50026" if pr.p < 0.05 else INK
                ax.errorbar(pr.r, y, xerr=[[pr.r - pr.lo], [pr.hi - pr.r]], fmt="D", color=col, ms=4.5, lw=1.3)
                ticks.append(y)
                labs.append(f"{v}: pooled (p = {pr.p:.2g})")
                y += 1.6
        ax.set_yticks(ticks, labs, fontsize=5.5)
        ax.invert_yaxis()
        ax.axvline(0, color=GRID, lw=0.8, zorder=0)
        ax.set_xlim(-0.8, 0.8)
        ax.set_xlabel("correlation with the axis within cases\n(duration, onset: adjusted for age at death)", fontsize=6)
        ax.set_title(AXIS_NAME[axn], fontsize=7, loc="left", fontweight="bold")
        ax.tick_params(labelsize=5.5, length=2)
        for sp_ in ["top", "right"]:
            ax.spines[sp_].set_visible(False)
    fig.subplots_adjust(left=0.24, right=0.98, top=0.92, bottom=0.16)
    save_figure(fig, FIG / "S5_axes_clinical")
    plt.close(fig)


if __name__ == "__main__":
    figures() if "--figs" in sys.argv else main()

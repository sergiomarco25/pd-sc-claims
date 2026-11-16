"""Side analysis (not part of the paper): do Lewy-body patients fall into reproducible molecular subgroups?

Usage: python -I scripts/subgroups/01_patient_subgroups.py
Patient profile: for every Lewy-body case (PD, PDD, DLB, LBD), programme and published-state scores per cell type
(results/biology/donor_scores.tsv.gz) and cell-type shares (donor_composition.tsv.gz), each expressed in standard
deviations from the controls of the same dataset x region (x cell type), then averaged over the regions of a donor
(the state is shared across regions of a patient). Cohorts analysed separately: those with >= 30 cases.
Per cohort: features with > 30% missing values are dropped and the rest imputed with 0 (the control mean); PCA on
standardised features; Gaussian mixtures (diagonal covariance) on the top 5 PCs for k = 1..5, compared by BIC
(k = 1 best = a continuum); k-means k = 2..4 with bootstrap stability (adjusted Rand index of 100 bootstrap
clusterings with the full-data clustering, on the shared donors).
Across cohorts (shared features): (1) replication of the main axes, absolute Tucker congruence of PC1-3 loadings
between cohorts, matched greedily; (2) replication of subgroups, k-means centroids of one cohort used to classify
another cohort's donors, compared with that cohort's own k-means (adjusted Rand index), and centroid profile
correlation. Clusters and axes are related to diagnosis, sex, age, Lewy stage and disease duration where known.
Writes results/subgroups/*.tsv and figures/subgroups/*.
"""
from pathlib import Path
import itertools
import sys
import warnings
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.cluster import KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score
from sklearn.mixture import GaussianMixture

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.plot.style import apply_style, save_figure, heatmap, CMAP_DIV, INK, MUTED, GRID  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402

warnings.simplefilter("ignore")

from pdrep import config as cfg  # noqa: E402

BASE = cfg.RESULTS
BIO, OUT = BASE / "biology", BASE / "subgroups"
FIG = ROOT / "figures" / "subgroups"
CASES = {"PD", "PDD", "DLB", "LBD"}
GLIA = ["Astrocytes", "Microglia", "Oligodendrocytes", "OPCs", "Endothelial"]
GEN = ["heat response", "chaperone folding", "unfolded protein response", "hypoxia", "cold-inducible (CIRBP, RBM3)",
       "inflammatory response", "interferon gamma", "oxidative phosphorylation"]
NEUR = [("heat response", c) for c in ["Excitatory neurons", "Inhibitory neurons"]] + \
       [("oxidative phosphorylation", c) for c in ["Excitatory neurons", "Inhibitory neurons"]] + \
       [("cold-inducible (CIRBP, RBM3)", c) for c in ["Excitatory neurons", "Inhibitory neurons"]]
STATES = [("disease-associated oligodendrocytes (Adams)", "Oligodendrocytes"),
          ("stress oligodendrocytes, Oligos5 (Martirosyan)", "Oligodendrocytes"),
          ("UPR astrocytes (Barba-Reyes)", "Astrocytes"), ("UPR astrocytes, Astrocytes3 (Martirosyan)", "Astrocytes"),
          ("C3/CD44 astrocytes, Astrocytes1 (Martirosyan)", "Astrocytes"), ("GFAP/MT astrocytes, Astrocytes4 (Martirosyan)", "Astrocytes"),
          ("AR2 astrocytes (Barba-Reyes)", "Astrocytes"), ("activated astrocytes, Ast-1 (Xu)", "Astrocytes"),
          ("UPR microglia (Barba-Reyes)", "Microglia"), ("GPNMB-high microglia (Smajić)", "Microglia"),
          ("HSP90AA1-high microglia (Smajić)", "Microglia"), ("activated microglia, Micr-1 (Xu)", "Microglia"),
          ("activated myeloid cells (Ma)", "Microglia"), ("DAM microglia (Barba-Reyes)", "Microglia")]
COMP = ["Microglia", "Astrocytes", "Oligodendrocytes", "OPCs", "Endothelial", "Excitatory neurons", "Inhibitory neurons",
        "neurons_total"]
RNG = np.random.default_rng(0)


def features():
    s = pd.read_csv(BIO / "donor_scores.tsv.gz", sep="\t", low_memory=False)
    want = {(p, c) for p in GEN for c in GLIA} | set(NEUR) | set(STATES)
    s = s[[(p, c) in want for p, c in zip(s.programme, s.cell_class)]].copy()

    def zc(d):
        r = d[d.disease.eq("control")].score
        return (d.score - r.mean()) / r.std(ddof=1) if len(r) >= 4 else d.score * np.nan
    s["z"] = s.groupby(["dataset", "region", "cell_class", "programme"], group_keys=False).apply(zc)
    s["feature"] = s.programme.str.replace(r" \(.*\)", "", regex=True) + " · " + s.cell_class
    meta = s.drop_duplicates(["dataset", "donor"])[["dataset", "donor", "disease", "sex", "age", "stage"]]
    sc = s[s.disease.isin(CASES)].groupby(["dataset", "donor", "feature"]).z.mean().unstack()
    c = pd.read_csv(BIO / "donor_composition.tsv.gz", sep="\t")
    c[["donor", "region"]] = c.unit.str.split("|", n=1, expand=True)
    c = c[c.cell_class.isin(COMP)].merge(meta[["dataset", "donor", "disease"]], on=["dataset", "donor"])
    c["logit"] = np.log((c.share + 1e-4) / (1 - c.share + 1e-4))

    def zc2(d):
        r = d[d.disease.eq("control")].logit
        return (d.logit - r.mean()) / r.std(ddof=1) if len(r) >= 4 else d.logit * np.nan
    c["z"] = c.groupby(["dataset", "region", "cell_class"], group_keys=False).apply(zc2)
    cc = c[c.disease.isin(CASES)].groupby(["dataset", "donor", "cell_class"]).z.mean().unstack()
    cc.columns = ["share · " + k for k in cc.columns]
    X = sc.join(cc, how="left")
    return X, meta.set_index(["dataset", "donor"])


def cohort_matrix(X, ds, max_missing=0.3):
    m = X.loc[ds]
    m = m.loc[:, m.isna().mean() <= max_missing]
    return m.fillna(0)


def analyse_cohort(m):
    Z = (m - m.mean()) / m.std(ddof=1).replace(0, 1)
    pca = PCA(n_components=min(10, Z.shape[1], Z.shape[0] - 1)).fit(Z)
    pcs = pca.transform(Z)[:, :5]
    bic = []
    for k in range(1, 6):
        g = GaussianMixture(k, covariance_type="diag", n_init=10, random_state=0).fit(pcs)
        bic.append(g.bic(pcs))
    stab = {}
    for k in (2, 3, 4):
        full = KMeans(k, n_init=20, random_state=0).fit_predict(Z)
        aris = []
        for _ in range(100):
            idx = RNG.choice(len(Z), len(Z), replace=True)
            u = np.unique(idx)
            lab = KMeans(k, n_init=5, random_state=int(RNG.integers(1e9))).fit(Z.values[idx]).predict(Z.values[u])
            aris.append(adjusted_rand_score(full[u], lab))
        stab[k] = (float(np.median(aris)), full)
    return Z, pca, pcs, bic, stab


def tucker(a, b):
    return float(np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b)))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    X, meta = features()
    X.to_csv(OUT / "patient_features.tsv.gz", sep="\t")
    n = X.groupby(level=0).size()
    cohorts = list(n[n >= 30].index)
    print("cases per cohort:", n.to_dict(), "\ncohorts analysed:", cohorts, flush=True)
    res, rows = {}, []
    for ds in cohorts:
        m = cohort_matrix(X, ds)
        Z, pca, pcs, bic, stab = analyse_cohort(m)
        res[ds] = (Z, pca, pcs, stab)
        best = int(np.argmin(bic)) + 1
        rows.append({"cohort": ds, "cases": len(Z), "features": Z.shape[1],
                     "pc1_var": pca.explained_variance_ratio_[0], "pc2_var": pca.explained_variance_ratio_[1],
                     "pc3_var": pca.explained_variance_ratio_[2], "best_k_bic": best,
                     **{f"bic_k{k + 1}": b for k, b in enumerate(bic)},
                     **{f"stability_k{k}": v[0] for k, v in stab.items()}})
    summ = pd.DataFrame(rows)
    summ.to_csv(OUT / "cohort_summary.tsv", sep="\t", index=False)
    print(summ.round(3).to_string(index=False), flush=True)
    # shared features across analysed cohorts
    shared = sorted(set.intersection(*[set(res[d][0].columns) for d in cohorts]))
    print("shared features:", len(shared), flush=True)
    load = {}
    for ds in cohorts:
        Z = res[ds][0][shared]
        p = PCA(n_components=3).fit(Z)
        load[ds] = pd.DataFrame(p.components_.T, index=shared, columns=["PC1", "PC2", "PC3"])
    axes = []
    for a, b in itertools.combinations(cohorts, 2):
        used = set()
        for i in range(3):
            best = max(((j, abs(tucker(load[a].iloc[:, i], load[b].iloc[:, j]))) for j in range(3) if j not in used),
                       key=lambda x: x[1])
            used.add(best[0])
            axes.append({"cohort_a": a, "cohort_b": b, "pc_a": i + 1, "pc_b": best[0] + 1, "congruence": best[1]})
    axes = pd.DataFrame(axes)
    axes.to_csv(OUT / "axis_replication.tsv", sep="\t", index=False)
    print(axes.round(2).to_string(index=False), flush=True)
    # subgroup replication via k-means centroids on shared features
    rep = []
    for k in (2, 3):
        cents = {}
        labs = {}
        for ds in cohorts:
            Z = res[ds][0][shared]
            km = KMeans(k, n_init=50, random_state=0).fit(Z)
            cents[ds], labs[ds] = km.cluster_centers_, km.labels_
        for a, b in itertools.permutations(cohorts, 2):
            Zb = res[b][0][shared].values
            pred = np.argmin(((Zb[:, None, :] - cents[a][None]) ** 2).sum(-1), axis=1)
            r = np.array([[stats.pearsonr(cents[a][i], cents[b][j])[0] for j in range(k)] for i in range(k)])
            rep.append({"k": k, "train": a, "test": b, "ari": adjusted_rand_score(labs[b], pred),
                        "best_centroid_r": float(np.mean(r.max(1)))})
    rep = pd.DataFrame(rep)
    rep.to_csv(OUT / "subgroup_replication.tsv", sep="\t", index=False)
    print(rep.groupby("k")[["ari", "best_centroid_r"]].describe().round(2).to_string(), flush=True)
    # relate axes to clinical variables
    clin = []
    dur = pd.read_csv(BASE / "survival" / "durations_used.tsv", sep="\t", dtype={"donor": str}).set_index(["dataset", "donor"])
    for ds in cohorts:
        Z, pca, pcs, stab = res[ds]
        mm = meta.loc[ds].reindex(Z.index)
        d = dur.duration.reindex([(ds, x) for x in Z.index]).values
        for i in range(3):
            pc = pcs[:, i]
            for var, v in [("age", pd.to_numeric(mm.age, errors="coerce").values),
                           ("Lewy stage", pd.to_numeric(mm.stage, errors="coerce").values), ("duration", d)]:
                ok = np.isfinite(v)
                if ok.sum() >= 15:
                    r, p = stats.spearmanr(pc[ok], v[ok])
                    clin.append({"cohort": ds, "pc": i + 1, "variable": var, "r": r, "p": p, "n": int(ok.sum())})
            for var, v in [("male", mm.sex.astype(str).str.upper().str[0].eq("M").values),
                           ("dementia (PDD/DLB)", mm.disease.isin(["PDD", "DLB"]).values)]:
                if 5 <= v.sum() <= len(v) - 5:
                    clin.append({"cohort": ds, "pc": i + 1, "variable": var, "r": np.nan,
                                 "p": stats.mannwhitneyu(pc[v], pc[~v]).pvalue, "n": len(v),
                                 "diff": pc[v].mean() - pc[~v].mean()})
    clin = pd.DataFrame(clin)
    clin.to_csv(OUT / "axes_vs_clinical.tsv", sep="\t", index=False)
    print(clin[clin.p < 0.05].round(3).to_string(index=False), flush=True)
    figures(res, cohorts, summ, load, shared)


def figures(res, cohorts, summ, load, shared):
    apply_style(font__size=7, axes__linewidth=0.6, axes__edgecolor=INK)
    fig, axs = plt.subplots(1, len(cohorts) + 1, figsize=(7.2, 2.4), gridspec_kw={"wspace": 0.45})
    for ax, ds in zip(axs, cohorts):
        pcs = res[ds][2]
        ax.scatter(pcs[:, 0], pcs[:, 1], s=5, color="#d4a400", edgecolor=INK, lw=0.2)
        v = res[ds][1].explained_variance_ratio_
        ax.set_xlabel(f"PC1 ({v[0]:.0%})", fontsize=6)
        ax.set_ylabel(f"PC2 ({v[1]:.0%})", fontsize=6)
        ax.set_title(f"{ds} (n = {len(pcs)})", fontsize=6.5, loc="left")
        ax.tick_params(labelsize=5.5)
    ax = axs[-1]
    for _, r in summ.iterrows():
        b = np.array([r[f"bic_k{k}"] for k in range(1, 6)])
        ax.plot(range(1, 6), b - b.min(), marker="o", ms=3, lw=0.8, label=r.cohort)
    ax.set_xlabel("mixture components k", fontsize=6)
    ax.set_ylabel("ΔBIC (0 = best)", fontsize=6)
    ax.legend(frameon=False, fontsize=5)
    ax.tick_params(labelsize=5.5)
    save_figure(fig, FIG / "S1_patient_space")
    plt.close(fig)
    # loadings of PC1-PC3 per cohort (shared features), signs aligned to the first cohort
    ref = load[cohorts[0]]
    mats, cols = [], []
    for ds in cohorts:
        L = load[ds].copy()
        for c in L.columns:
            if np.dot(L[c], ref[c]) < 0:
                L[c] = -L[c]
        mats.append(L)
        cols += [f"{ds}\n{c}" for c in L.columns]
    M = pd.concat(mats, axis=1)
    order = M.iloc[:, 0].sort_values().index
    fig, ax = plt.subplots(figsize=(7.2, 0.12 * len(order) + 1.2))
    qm = heatmap(ax, M.loc[order].values, cmap=CMAP_DIV, vmin=-0.35, vmax=0.35)
    ax.set_yticks(range(len(order)), order, fontsize=5)
    ax.set_xticks(range(M.shape[1]), cols, fontsize=5, rotation=90)
    cb = fig.colorbar(qm, ax=ax, fraction=0.02, pad=0.01)
    cb.set_label("PCA loading", fontsize=6)
    save_figure(fig, FIG / "S2_axis_loadings")
    plt.close(fig)


if __name__ == "__main__":
    main()

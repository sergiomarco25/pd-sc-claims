"""Side analysis: interpret the patient axes that replicate across cohorts (scripts/subgroups/01_patient_subgroups.py).

Usage: python -I scripts/subgroups/02_axes_interpretation.py
PCA per cohort on the shared features (signs aligned to Menon 2026); mean loading of each feature over the four cohorts
for PC1 and PC2; correlation of each cohort's PC1/PC2 with RNA quality (log median genes per nucleus over the donor's
units, results/biology/qc_units.tsv.gz) and post-mortem interval. Writes results/subgroups/axis_loadings.tsv and
axes_vs_quality.tsv.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats
from sklearn.decomposition import PCA

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

R = cfg.RESULTS
COHORTS = ["barbareyes2026", "lee2026", "menon2026", "prashant2024"]


def consensus_loadings(X, cohorts, ref="menon2026"):
    """Case profiles X (dataset, donor) x feature -> PC1/PC2 loadings per cohort on the features all cohorts share
    (each cohort: features > 30% missing dropped, rest 0, standardised), signs aligned to `ref`, plus their mean;
    and each cohort's PC scores."""
    Zs = {}
    for ds in cohorts:
        m = X.loc[ds]
        m = m.loc[:, m.isna().mean() <= 0.3].fillna(0)
        Zs[ds] = (m - m.mean()) / m.std(ddof=1).replace(0, 1)
    shared = sorted(set.intersection(*[set(z.columns) for z in Zs.values()]))
    L, P = {}, {}
    for ds in cohorts:
        p = PCA(3).fit(Zs[ds][shared])
        L[ds] = pd.DataFrame(p.components_.T, index=shared, columns=["PC1", "PC2", "PC3"])
        P[ds] = pd.DataFrame(p.transform(Zs[ds][shared]), index=Zs[ds].index, columns=["PC1", "PC2", "PC3"])
    for ds in cohorts:
        for c in ["PC1", "PC2"]:
            if np.dot(L[ds][c], L[ref][c]) < 0:
                L[ds][c] *= -1
                P[ds][c] *= -1
    load = pd.concat({ds: L[ds][["PC1", "PC2"]] for ds in cohorts}, axis=1)
    load[("mean", "PC1")] = load.xs("PC1", axis=1, level=1).mean(1)
    load[("mean", "PC2")] = load.xs("PC2", axis=1, level=1).mean(1)
    return load, P


def main():
    X = pd.read_csv(R / "subgroups" / "patient_features.tsv.gz", sep="\t", index_col=[0, 1])
    load, P = consensus_loadings(X, COHORTS)
    load.to_csv(R / "subgroups" / "axis_loadings.tsv", sep="\t")
    q = pd.read_csv(R / "biology" / "qc_units.tsv.gz", sep="\t")
    q = q[q.cell_class == "ALL"].assign(donor=lambda d: d.unit.str.split("|").str[0])
    q = q.groupby(["dataset", "donor"]).median_genes.apply(lambda v: np.log(v).mean())
    pm = pd.read_csv(R / "biology" / "donor_scores.tsv.gz", sep="\t", usecols=["dataset", "donor", "pmi"],
                     low_memory=False).drop_duplicates(["dataset", "donor"]).set_index(["dataset", "donor"]).pmi
    rows = []
    for ds in COHORTS:
        for var, v in [("log genes per nucleus", q), ("post-mortem interval", pm)]:
            x = pd.to_numeric(v.reindex([(ds, d) for d in P[ds].index]), errors="coerce").values
            ok = np.isfinite(x)
            for c in ["PC1", "PC2"]:
                if ok.sum() > 10:
                    r, p = stats.spearmanr(P[ds][c].values[ok], x[ok])
                    rows.append({"cohort": ds, "axis": c, "variable": var, "r": r, "p": p, "n": int(ok.sum())})
    pd.DataFrame(rows).to_csv(R / "subgroups" / "axes_vs_quality.tsv", sep="\t", index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False))


if __name__ == "__main__":
    main()

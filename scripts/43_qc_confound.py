"""Is the stress score an RNA-quality artefact? Relation to per-donor nucleus quality, and effects adjusted for it.

Usage: python -I scripts/43_qc_confound.py
For every harmonised dataset, per donor x region x cell class: number of nuclei, median genes and median counts per
nucleus (from obs). Joined to the programme scores (results/biology/donor_scores.tsv.gz):
  1. within each dataset x region x class, Spearman correlation of the score with log median genes per nucleus, in
     cases and in controls, pooled across datasets (DerSimonian-Laird on Fisher z);
  2. case-control Hedges' g before and after regressing the score on log median genes and log nuclei within the
     dataset x region x class (residuals), pooled across datasets;
  3. case-control difference in log median genes itself.
Writes results/biology/qc_units.tsv.gz and qc_confound.tsv.
"""
from pathlib import Path
import sys
import warnings
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.meta import fixed, dersimonian_laird  # noqa: E402
from pdrep.harmonise import OUT as H  # noqa: E402

warnings.simplefilter("ignore")

from pdrep import config as cfg  # noqa: E402

B = cfg.RESULTS / "biology"
CASES = {"PD", "PDD", "DLB", "LBD"}
PROGS = [("heat response", "Oligodendrocytes"), ("heat response", "Astrocytes"), ("heat response", "Microglia"),
         ("chaperone folding", "Oligodendrocytes"), ("cold-inducible (CIRBP, RBM3)", "Oligodendrocytes"),
         ("cold-inducible (CIRBP, RBM3)", "ALL"), ("heat response", "ALL"), ("hypoxia", "Microglia"),
         ("UPR microglia (Barba-Reyes)", "Microglia"), ("UPR astrocytes, Astrocytes3 (Martirosyan)", "Astrocytes"),
         ("disease-associated oligodendrocytes (Adams)", "Oligodendrocytes")]


def qc_units():
    rows = []
    for f in sorted(H.glob("*.h5ad")):
        o = ad.read_h5ad(f, backed="r").obs
        if "enrichment" in o:
            o = o[~o.enrichment.astype(str).str.lower().eq("sorted")]
        g = o.assign(unit=o.donor.astype(str) + "|" + o.region.astype(str))
        for cls_, d in list(g.groupby("cell_class", observed=True)) + [("ALL", g)]:
            q = d.groupby("unit", observed=True).agg(nuclei=("n_genes_by_counts", "size"),
                                                      median_genes=("n_genes_by_counts", "median"),
                                                      median_counts=("total_counts", "median")).reset_index()
            rows.append(q.assign(dataset=f.stem, cell_class=str(cls_)))
        print(f.stem, flush=True)
    q = pd.concat(rows, ignore_index=True)
    q.to_csv(B / "qc_units.tsv.gz", sep="\t", index=False)
    return q


def hedges(a, r):
    if len(a) < 3 or len(r) < 3:
        return np.nan, np.nan
    sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(r) - 1) * r.var(ddof=1)) / (len(a) + len(r) - 2))
    if not sp > 0:
        return np.nan, np.nan
    g = (a.mean() - r.mean()) / sp
    return g, np.sqrt(1 / len(a) + 1 / len(r) + g ** 2 / (2 * (len(a) + len(r))))


def pool(rows, key):
    out = []
    d = pd.DataFrame(rows).dropna(subset=["est"])
    for k, x in d.groupby(key):
        per = x.groupby("dataset").apply(lambda y: pd.Series(fixed(y.est, y.se), index=["est", "se"]))
        r = dersimonian_laird(per.est, per.se)
        out.append({**dict(zip(key, k)), **r})
    return pd.DataFrame(out)


def main():
    q = qc_units() if not (B / "qc_units.tsv.gz").exists() else pd.read_csv(B / "qc_units.tsv.gz", sep="\t")
    s = pd.read_csv(B / "donor_scores.tsv.gz", sep="\t", low_memory=False)
    s["key"] = s.donor.astype(str) + "|" + s.region.astype(str)
    m = s.merge(q, left_on=["dataset", "cell_class", "key"], right_on=["dataset", "cell_class", "unit"])
    m["lg"], m["ln"] = np.log(m.median_genes), np.log(m.nuclei)
    rows = []
    for (p, c) in PROGS:
        d = m[(m.programme == p) & (m.cell_class == c)]
        for (ds, reg), x in d.groupby(["dataset", "region"]):
            x = x.dropna(subset=["score", "lg"])
            cas, con = x[x.disease.isin(CASES)], x[x.disease.eq("control")]
            for grp, y in [("cases", cas), ("controls", con)]:
                if len(y) >= 8:
                    r = stats.spearmanr(y.score, y.lg)[0]
                    rows.append({"programme": p, "cell_class": c, "test": f"score ~ genes per nucleus ({grp})",
                                 "dataset": ds, "est": np.arctanh(np.clip(r, -0.999, 0.999)), "se": 1 / np.sqrt(len(y) - 3)})
            y = pd.concat([cas, con])
            if len(cas) >= 3 and len(con) >= 3:
                g, se = hedges(cas.score.values, con.score.values)
                rows.append({"programme": p, "cell_class": c, "test": "cases vs controls", "dataset": ds, "est": g, "se": se})
                X = np.c_[np.ones(len(y)), y.lg, y.ln]
                beta = np.linalg.lstsq(X, y.score.values, rcond=None)[0]
                res = y.score.values - X @ beta
                isc = y.disease.isin(CASES).values
                g, se = hedges(res[isc], res[~isc])
                rows.append({"programme": p, "cell_class": c, "test": "cases vs controls, adjusted for nucleus quality",
                             "dataset": ds, "est": g, "se": se})
                g, se = hedges(cas.lg.values, con.lg.values)
                rows.append({"programme": p, "cell_class": c, "test": "genes per nucleus, cases vs controls",
                             "dataset": ds, "est": g, "se": se})
    out = pool(rows, ["programme", "cell_class", "test"])
    corr = out.test.str.startswith("score ~")
    out.loc[corr, "r"] = np.tanh(out.loc[corr, "est"])
    out.to_csv(B / "qc_confound.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(out[["programme", "cell_class", "test", "est", "r", "p", "k"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()

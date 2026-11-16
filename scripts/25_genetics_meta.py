"""Which cell types carry PD genetic risk consistently? Meta-analysis of the per-dataset MAGMA tests.

Usage: python -I scripts/25_genetics_meta.py
Inputs: results/gwas/sets/<dataset>.gsa.out (MAGMA competitive tests of top-10% cell-class-specific genes and of PD-up
genes, Nalls 2019 GWAS) and <dataset>.familial.tsv (familial PD gene overlap). Per cell class and set type, one-sided
z (= Phi^-1(1 - p)) are combined across datasets by Stouffer; also the share of datasets with p < 0.05. Familial
overlap tests (discrete, p = 1 for empty overlaps) are combined by a binomial test of that share against 5%.
Writes results/biology/genetics_meta.tsv and genetics_per_dataset.tsv.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

BASE = cfg.RESULTS
SETS, OUT = BASE / "gwas" / "sets", BASE / "biology"


def pool(per):
    """Per-dataset tests -> per set type x cell class: Stouffer z of MAGMA tests, binomial test for familial overlap."""
    g = per.groupby(["set_type", "cell_class"])
    meta = pd.DataFrame({"k": g.size(), "z": g.z.sum() / np.sqrt(g.size()), "share_p05": g.p.apply(lambda p: (p < 0.05).mean()),
                         "datasets_p05": g.apply(lambda d: ";".join(sorted(d.dataset[d.p < 0.05])))}).reset_index()
    meta["p"] = stats.norm.sf(meta.z)
    # familial overlap tests are discrete (p = 1 whenever a set holds no familial gene), so Stouffer is not suitable:
    # combine them by the number of datasets with p < 0.05 against the 5% expected by chance (binomial)
    fam = meta.set_type.eq("familial PD genes")
    nsig = (meta.share_p05 * meta.k).round().astype(int)
    meta.loc[fam, "p"] = stats.binom.sf(nsig[fam] - 1, meta.k[fam], 0.05)
    meta.loc[fam, "z"] = stats.norm.isf(meta.loc[fam, "p"].clip(1e-300, 1 - 1e-12))
    return meta


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for f in sorted(SETS.glob("*.gsa.out")):
        ds = f.name[: -len(".gsa.out")]
        t = pd.read_csv(f, sep=r"\s+", comment="#")
        name = (t["FULL_NAME"] if "FULL_NAME" in t else t["VARIABLE"]).astype(str)
        for n, p, b, k in zip(name, t.P, t.BETA, t.NGENES):
            kind, *rest = n.split("::")
            rows.append({"dataset": ds, "set_type": "specificity" if kind == "spec" else "PD-up genes (" + rest[0] + ")",
                         "cell_class": rest[-1].replace("_", " ").replace("Immune other", "Immune (other)")
                         .replace("Neurons other", "Neurons (other)"), "p": p, "beta": b, "n_genes": k, "test": "MAGMA PD GWAS"})
        fam = SETS / f"{ds}.familial.tsv"
        if fam.exists():
            for r in pd.read_csv(fam, sep="\t").itertuples():
                rows.append({"dataset": ds, "set_type": "familial PD genes", "cell_class": r.cell_class, "p": r.p,
                             "beta": r.familial_in_set, "n_genes": r.n_genes, "test": "hypergeometric"})
    per = pd.DataFrame(rows)
    per["z"] = stats.norm.isf(per.p.clip(1e-300, 1 - 1e-12))
    per = per[per.cell_class != "ALL"]
    per.to_csv(OUT / "genetics_per_dataset.tsv", sep="\t", index=False)
    meta = pool(per)
    meta.to_csv(OUT / "genetics_meta.tsv", sep="\t", index=False)
    print(meta[meta.k >= 3].sort_values(["set_type", "p"]).round(4).to_string(index=False))


if __name__ == "__main__":
    main()

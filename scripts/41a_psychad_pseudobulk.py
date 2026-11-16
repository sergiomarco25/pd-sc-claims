"""Pseudobulk counts for every PsychAD donor (Lee et al. 2026; DLPFC; 1,042 donors), all diagnoses.

Usage: python -I scripts/41a_psychad_pseudobulk.py
The harmonised lee2026 dataset keeps only Lewy-body donors and neurotypical controls. To ask whether the PD stress
state is specific to Lewy-body disease, this script streams the full CELLxGENE file (4.1 million nuclei) in row chunks
and sums raw counts per donor x cell class (classes from the author subtypes, as in src/pdrep/build/lee2026.py),
without loading the matrix. It also records per donor x class the number of nuclei and the median genes and counts per
nucleus (RNA-quality covariates).
Writes results/psychad/pb__<class>.csv.gz (genes x donors) and results/psychad/donors.tsv (diagnosis flags, groups).
Diagnostic groups (exclusive, by the release flags): control (neurotypical); Lewy-body disease without AD; Lewy-body
disease with AD; AD (no Lewy-body disease); vascular dementia only; FTD; tauopathy; schizophrenia without dementia;
bipolar disorder without dementia; other dementia.
"""
from pathlib import Path
import sys
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.build._io import read_df  # noqa: E402
from pdrep.build.lee2026 import PATH, _cls, _age  # noqa: E402

from pdrep import config as cfg  # noqa: E402

OUT = cfg.RESULTS / "psychad"
CLASSES = ["Oligodendrocytes", "OPCs", "Astrocytes", "Microglia", "Excitatory neurons", "Inhibitory neurons", "Endothelial"]


def group(r):
    lewy = r.DLBD_status == "Yes" or r.Parkinson_disease == "Yes"
    ad = r.AD_status == "Yes"
    terms = set(str(r.disease).split(" || "))
    if r.disease == "normal":
        return "control"
    if lewy:
        return "Lewy-body disease with AD" if ad else "Lewy-body disease"
    if ad:
        return "Alzheimer's disease"
    if r.FTD_status == "Yes":
        return "frontotemporal dementia"
    if r.Tauopathy_status == "Yes":
        return "tauopathy"
    if r.Vascular_status == "Yes":
        return "vascular dementia"
    if "schizophrenia" in terms:
        return "schizophrenia" if "dementia" not in terms else "schizophrenia with dementia"
    if r.Bipolar_Disorder == "Yes" or "bipolar disorder" in terms:
        return "bipolar disorder" if "dementia" not in terms else "other dementia"
    return "other dementia" if "dementia" in terms else "other"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    o = read_df(PATH, "obs")
    var = read_df(PATH, "var")
    genes = var["feature_name"].astype(str).values
    cls = o["subtype"].astype(str).map(_cls).values
    donor = o["donor_id"].astype(str).values
    donors = pd.Index(sorted(set(donor)))
    d = o.drop_duplicates("donor_id").set_index("donor_id")
    meta = pd.DataFrame({"donor": donors})
    for c in ["disease", "AD_status", "DLBD_status", "FTD_status", "Tauopathy_status", "Vascular_status",
              "Schizophrenia", "Bipolar_Disorder", "Parkinson_disease", "sex"]:
        meta[c] = d.loc[donors, c].astype(str).values
    meta["age"] = d.loc[donors, "development_stage"].astype(str).map(_age).values
    meta["group"] = [group(r) for r in meta.itertuples()]
    meta.to_csv(OUT / "donors.tsv", sep="\t", index=False)
    print(meta.group.value_counts().to_string(), flush=True)
    di = donors.get_indexer(donor)
    sums = {c: np.zeros((len(genes), len(donors)), np.float64) for c in CLASSES}
    qc = []
    with h5py.File(PATH, "r") as h:
        g = h["X"]
        n_rows, n_cols = g.attrs["shape"]
        indptr = g["indptr"][:]
        for s in range(0, n_rows, 200_000):
            e = min(s + 200_000, n_rows)
            a, b = indptr[s], indptr[e]
            m = sp.csr_matrix((g["data"][a:b], g["indices"][a:b], indptr[s:e + 1] - a), shape=(e - s, n_cols))
            for c in CLASSES:
                rows = np.where(cls[s:e] == c)[0]
                if not len(rows):
                    continue
                ind = sp.csr_matrix((np.ones(len(rows)), (di[s:e][rows], np.arange(len(rows)))),
                                    shape=(len(donors), len(rows)))
                sums[c] += (ind @ m[rows]).T.toarray()
            print(e, flush=True)
    obs_qc = pd.DataFrame({"donor": donor, "cell_class": cls, "n_genes": o["n_genes"].values,
                           "n_counts": o["n_counts"].values})
    q = obs_qc[obs_qc.cell_class.isin(CLASSES)].groupby(["donor", "cell_class"]).agg(
        nuclei=("n_genes", "size"), median_genes=("n_genes", "median"), median_counts=("n_counts", "median"))
    q.reset_index().to_csv(OUT / "qc.tsv", sep="\t", index=False)
    for c in CLASSES:
        df = pd.DataFrame(sums[c], index=genes, columns=donors)
        df = df.groupby(level=0).sum()
        df = df.loc[:, df.sum(0) > 0]
        df.round().astype(np.int64).to_csv(OUT / f"pb__{c.replace(' ', '_')}.csv.gz")
        print(c, df.shape, flush=True)


if __name__ == "__main__":
    main()

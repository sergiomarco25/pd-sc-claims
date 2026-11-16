"""lee2023 — Lee et al. 2023 Sci Adv, GSE148434 (10x v3 snRNA-seq, substantia nigra, UCSD ADRC; in-house donors only).

Sources (raw/lee2023/):
  GSM*_<donor>_snRNA.molecule_info.h5   Cell Ranger molecule_info (file_version 3), one per donor (12: NOSN01-06,
                                         PDSN01/02/04/06/07/09). No count matrices are deposited.
  supplements/lee2023/tables_s1_to_s10/abo2467_Table_S1.xlsx   donor metadata (diagnosis, APOE, age, sex, PMI,
                                         "Braak combo", MMSE, DRS, snRNA-seq pool number).
Counts: per-cell UMI counts are rebuilt from molecule_info: every molecule entry (one per unique barcode x UMI x
feature after Cell Ranger's UMI correction/dedup) with a barcode in barcode_info/pass_filter adds 1 to
X[barcode, feature]. Pooled libraries (pools 1-4: NOSN01/02/04, NOSN05/06, PDSN01/02, PDSN04/06/07) were deposited
already demultiplexed (souporcell): their molecule_info files contain only the donor's pass_filter barcodes
(checked: no barcode shared between donors of the same pool). The two unpooled libraries (NOSN03, PDSN09) contain all
barcodes and are restricted to pass_filter. Then uniform QC: >= 200 genes, < 10% mito UMIs.
The 7 external controls the paper re-processed from GSE140231 (Agarwal 2020) are not included here.
disease: PD / control. braak_lb NaN: Table S1 gives only a "Braak combo" value (e.g. 4.2), kept as `braak_combo`.
No public author cell labels -> marker clustering in finalise().
"""
import re
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from pdrep.harmonise import RAW

AUTHOR_TO_CLASS = {}
D = RAW / "lee2023"
S1 = RAW / "supplements" / "lee2023" / "tables_s1_to_s10" / "abo2467_Table_S1.xlsx"


def table_s1():
    v = pd.read_excel(S1, header=None)
    t = v.iloc[3:23, :14].copy()
    t.columns = ["source", "donor", "sample_id", "diagnosis", "region", "apoe", "death", "age", "sex", "pmi",
                 "braak_combo", "mmse", "drs", "snrna_pool"]
    t = t[t.donor.notna()].set_index("donor")
    return t


def read_molinfo(path):
    with h5py.File(path, "r") as f:
        bcs = f["barcodes"][:].astype(str)
        pf = np.unique(f["barcode_info/pass_filter"][:, 0])
        bi = f["barcode_idx"][:]
        fi = f["feature_idx"][:]
        names = f["features/name"][:].astype(str)
        ids = f["features/id"][:].astype(str)
        ftype = f["features/feature_type"][:].astype(str)
    keep = np.isin(bi, pf)
    row = np.searchsorted(pf, bi[keep])
    X = sp.csr_matrix((np.ones(keep.sum(), dtype=np.float32), (row, fi[keep])), shape=(len(pf), len(names)))
    X.sum_duplicates()
    g = ftype == "Gene Expression"
    a = ad.AnnData(X[:, g], obs=pd.DataFrame(index=bcs[pf] + "-1"),
                   var=pd.DataFrame({"gene_ids": ids[g]}, index=names[g]))
    a.var_names_make_unique()
    return a, int((~keep).sum())


def qc(a, min_genes=200, max_mito=10.0):
    X = a.X.tocsr()
    ng = np.diff(X.indptr)
    tot = np.asarray(X.sum(1)).ravel()
    mt = np.asarray(X[:, a.var_names.str.startswith("MT-")].sum(1)).ravel()
    keep = (ng >= min_genes) & (100 * mt / np.maximum(tot, 1) < max_mito)
    return a[keep].copy()


def load():
    meta = table_s1()
    parts = []
    for f in sorted(D.glob("GSM*_snRNA.molecule_info.h5")):
        gsm, donor = re.match(r"(GSM\d+)_(\w+SN\d+)_snRNA", f.name).groups()
        a, n_out = read_molinfo(f)
        n0 = a.n_obs
        a = qc(a)
        print(f"{donor}: pass_filter {n0} (molecules outside pass_filter dropped: {n_out}), after QC {a.n_obs}, "
              f"median UMI {np.median(np.asarray(a.X.sum(1)).ravel()):.0f}", flush=True)
        m = meta.loc[donor]
        a.obs_names = donor + "_" + a.obs_names
        a.obs["sample"] = donor
        a.obs["gsm"] = gsm
        a.obs["donor"] = donor
        a.obs["region"] = "SN"
        a.obs["disease"] = {"Normal": "control", "Parkinson": "PD"}[m["diagnosis"]]
        a.obs["disease_raw"] = m["diagnosis"]
        a.obs["sex"] = {"Female": "F", "Male": "M"}[m["sex"]]
        a.obs["age"] = float(m["age"])
        a.obs["pmi"] = pd.to_numeric(m["pmi"], errors="coerce")
        a.obs["braak_lb"] = np.nan
        a.obs["braak_combo"] = str(m["braak_combo"])
        a.obs["mmse"] = pd.to_numeric(m["mmse"], errors="coerce")
        a.obs["drs"] = pd.to_numeric(m["drs"], errors="coerce")
        a.obs["apoe"] = str(m["apoe"])
        a.obs["snrna_pool"] = str(m["snrna_pool"])
        a.obs["author_cell_type"] = "unknown"
        parts.append(a)
    return ad.concat(parts, join="inner", merge="same")

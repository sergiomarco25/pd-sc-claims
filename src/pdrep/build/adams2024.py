"""adams2024 — Adams et al. 2024 Nat Aging, GSE193688 (10x Multiome, midbrain / substantia nigra).

Sources (raw/adams2024/):
  GSM*_<sample>_filtered_feature_bc_matrix.h5   Cell Ranger ARC 1.0.1 filtered matrices (34 libraries, one per donor);
                                                only the "Gene Expression" features are read (ATAC peaks dropped).
  GSE193688_family.soft.gz                      group (Young/Aged/PD), diagnosis, age, sex, PMI, RIN, repository.
  supplements/adams2024/..._MOESM3_ESM.xlsx     Supp. Table 1: same metadata + 3 donors flagged
                                                "Poor RNA Quality, Omitted from downstream analysis" (YC5442, C4634, P6320).
Counts: integer UMI counts straight from the filtered matrices (Cell Ranger ARC cell calling).
QC here (uniform across loaders): >= 200 genes and < 10% mitochondrial UMIs per nucleus; the three donors the
authors omitted are dropped. No author cell labels are public -> cells labelled by marker clustering in finalise().
Region: "midbrain" (blocks dissected around the SN; GEO source "substantia nigra"); region_raw keeps the source text.
age_group: young / aged / PD (both young and aged controls have disease = control).
"""
import gzip
import re
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad
from pdrep.harmonise import RAW

AUTHOR_TO_CLASS = {}
D = RAW / "adams2024"
OMITTED = {"YC5442", "C4634", "P6320"}  # Supp. Table 1: poor RNA quality, omitted by the authors


def soft_samples(path):
    rows, cur = [], None
    with gzip.open(path, "rt", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                cur = {"gsm": line.split("=", 1)[1].strip()}
                rows.append(cur)
            elif cur is not None and line.startswith("!Sample_title"):
                cur["title"] = line.split("=", 1)[1].strip()
            elif cur is not None and line.startswith("!Sample_source_name_ch1"):
                cur["source"] = line.split("=", 1)[1].strip()
            elif cur is not None and line.startswith("!Sample_characteristics_ch1"):
                k, _, v = line.split("=", 1)[1].strip().partition(":")
                cur[k.strip().lower()] = v.strip()
    return pd.DataFrame(rows)


def qc(a, min_genes=200, max_mito=10.0):
    a.var["mt"] = a.var_names.str.startswith("MT-")
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    keep = (a.obs.n_genes_by_counts >= min_genes) & (a.obs.pct_counts_mt < max_mito)
    a = a[keep].copy()
    a.obs = a.obs.drop(columns=[c for c in a.obs.columns if c.endswith("counts") or "counts_mt" in c or c == "n_genes_by_counts"])
    a.var = a.var.drop(columns=[c for c in a.var.columns if c not in ("gene_ids",)])
    return a


def load():
    meta = soft_samples(D / "GSE193688_family.soft.gz").set_index("title")
    parts, n_raw = [], {}
    for f in sorted(D.glob("GSM*_filtered_feature_bc_matrix.h5")):
        s = re.match(r"GSM\d+_(.+)_filtered_feature_bc_matrix\.h5", f.name).group(1)
        if s in OMITTED:
            continue
        a = sc.read_10x_h5(f, gex_only=True)
        a.var_names_make_unique()
        a.var = a.var[["gene_ids"]]
        n_raw[s] = a.n_obs
        a.obs_names = s + "_" + a.obs_names
        a.obs["sample"] = s
        parts.append(qc(a))
    a = ad.concat(parts, join="inner", merge="same")
    print("nuclei per sample before/after QC:", {s: (n_raw[s], int((a.obs["sample"] == s).sum())) for s in n_raw})
    m = meta.loc[a.obs["sample"].values]
    a.obs["donor"] = a.obs["sample"].values
    a.obs["region"] = "midbrain"
    a.obs["region_raw"] = m["source"].values
    grp = m["group"].str.lower().values
    a.obs["age_group"] = grp
    a.obs["disease"] = np.where(grp == "pd", "PD", "control")
    a.obs["disease_raw"] = (m["group"] + ": " + m["diagnosis"]).values
    a.obs["sex"] = m["sex"].values
    a.obs["age"] = pd.to_numeric(m["age"], errors="coerce").values
    a.obs["pmi"] = pd.to_numeric(m["pmi (hours)"], errors="coerce").values
    a.obs["rin"] = pd.to_numeric(m["rin"], errors="coerce").values
    a.obs["brain_bank"] = m["repository"].values
    a.obs["braak_lb"] = np.nan
    a.obs["author_cell_type"] = "unknown"
    return a

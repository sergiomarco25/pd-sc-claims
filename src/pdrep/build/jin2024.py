"""jin2024 — Jin et al. 2024 Sci Adv, GSE235914 (10x 3' v3.1 snRNA-seq, superior temporal cortex, Mayo Clinic).

Sources (raw/jin2024/):
  GSE235914_RAW.tar -> extracted/GSM*_<n>-filtered-{matrix.mtx,features.tsv,barcodes.tsv}.gz  (only the
      Cell Ranger *filtered* matrices were extracted; the raw/unfiltered ones are ignored)
  GSE235914_family.soft.gz   npid (donor), type (Ctrl/LBD), pathdx (TLBD / DLBD / Dupl_DLBD / Tri_DLBD), sex, age,
      APOE, RIN, batch, sequencing saturation. Titles/characteristics of GSM7511853 (Tri) and GSM7511857 (Dup)
      were corrected by the submitters on 2026-03-04; the current soft file is used.
  supplements/jin2024/sciadv.adk3700_sm.pdf  Table S4: SNCA allele number, Braak (NFT) stage, Thal phase per donor
      (matched here by type + age, which is unique; ages agree with GEO).
Counts: integer UMI counts from the filtered matrices. QC: >= 200 genes, < 10% mito UMIs per nucleus.
No public author cell labels -> marker clustering in finalise().
disease: all LBD donors -> "LBD" (clinical PD/PDD/DLB subtype not given); pathdx kept in disease_raw.
genotype: SNCA_normal / SNCA_dup / SNCA_tri (from pathdx, consistent with Table S4 allele numbers 2/3/4).
braak_lb: NaN (Table S4 Braak stage is the NFT stage -> braak_nft).
"""
import gzip
import re
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad
import scanpy as sc
from pdrep.harmonise import RAW

AUTHOR_TO_CLASS = {}
D = RAW / "jin2024"
E = D / "extracted"
# Table S4 (supplementary PDF): (type, age) -> (SNCA alleles, Braak NFT, Thal)
TABLE_S4 = {("Ctrl", 68): (2, 1, 0), ("Ctrl", 67): (2, 1, 0), ("Ctrl", 59): (2, 0, 0), ("LBD", 68): (2, 2, 1),
            ("LBD", 55): (2, 3, 0), ("LBD", 66): (2, 3, 0), ("LBD", 63): (3, 1, 2), ("LBD", 48): (4, 1, 0)}


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
            elif cur is not None and line.startswith("!Sample_characteristics_ch1"):
                k, _, v = line.split("=", 1)[1].strip().partition(":")
                cur[k.strip().lower()] = v.strip()
    return pd.DataFrame(rows).set_index("gsm")


def qc(a, min_genes=200, max_mito=10.0):
    X = a.X.tocsr()
    ng = np.diff(X.indptr)
    tot = np.asarray(X.sum(1)).ravel()
    mt = np.asarray(X[:, a.var_names.str.startswith("MT-")].sum(1)).ravel()
    keep = (ng >= min_genes) & (100 * mt / np.maximum(tot, 1) < max_mito)
    return a[keep].copy()


def read_mtx(prefix):
    X = sp.csr_matrix(scipy.io.mmread(str(prefix) + "matrix.mtx.gz").T.tocsr())
    feat = pd.read_csv(str(prefix) + "features.tsv.gz", sep="\t", header=None)
    bc = pd.read_csv(str(prefix) + "barcodes.tsv.gz", sep="\t", header=None)[0].values
    feat = feat[feat[2] == "Gene Expression"] if feat.shape[1] > 2 else feat
    a = ad.AnnData(X[:, feat.index.values], obs=pd.DataFrame(index=bc),
                   var=pd.DataFrame({"gene_ids": feat[0].values}, index=feat[1].values))
    a.var_names_make_unique()
    return a


def load():
    meta = soft_samples(D / "GSE235914_family.soft.gz")
    parts = []
    for f in sorted(E.glob("GSM*-filtered-matrix.mtx.gz")):
        gsm, n = re.match(r"(GSM\d+)_(\d+)-filtered-matrix", f.name).groups()
        a = read_mtx(E / f"{gsm}_{n}-filtered-")
        n0 = a.n_obs
        a = qc(a)
        print(gsm, "nuclei before/after QC", n0, a.n_obs, flush=True)
        m = meta.loc[gsm]
        age = int(float(m["age"]))
        snca, braak_nft, thal = TABLE_S4[(m["type"], age)]
        a.obs_names = f"{gsm}_" + a.obs_names
        a.obs["sample"] = gsm
        a.obs["donor"] = m["npid"]
        a.obs["sample_title"] = m["title"]
        a.obs["region"] = "temporal"
        a.obs["disease"] = "LBD" if m["type"] == "LBD" else "control"
        a.obs["disease_raw"] = f'{m["type"]}: {m["pathdx"]}'
        a.obs["genotype"] = {2: "SNCA_normal", 3: "SNCA_dup", 4: "SNCA_tri"}[snca]
        a.obs["snca_alleles"] = snca
        a.obs["sex"] = m["sex"]
        a.obs["age"] = float(age)
        a.obs["pmi"] = np.nan
        a.obs["braak_lb"] = np.nan
        a.obs["braak_nft"] = float(braak_nft)
        a.obs["thal"] = float(thal)
        a.obs["apoe"] = m["apoe"]
        a.obs["rin"] = float(m["rin"])
        a.obs["batch"] = m["batch"]
        a.obs["author_cell_type"] = "unknown"
        parts.append(a)
    # sanity: GEO pathdx and Table S4 allele numbers must agree
    a = ad.concat(parts, join="inner", merge="same")
    dup = a.obs.drop_duplicates("donor")
    assert ((dup.disease_raw.str.contains("Dupl") == (dup.snca_alleles == 3)) &
            (dup.disease_raw.str.contains("Tri") == (dup.snca_alleles == 4))).all()
    return a

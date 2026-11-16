"""glock2026 — Glock et al. 2026, GSE295998 (10x snRNA-seq, amygdala; 15 late-onset PD, 8 controls, 24 libraries).

Sources (raw/glock2026/):
  GSM*_<lib>_cellbender_filtered.h5   CellBender v0.2 (default params) output of Cell Ranger v6 raw matrices
                                      (GRCh38, introns included), one per library (L1-8, K1-8, J1-8).
  GSE295998_UmapCoordinatesAndClusterAssignments.txt.gz  authors' post-QC nuclei (168,270, all scDblFinder
                                      singlets) with Library/Subject/Condition/Pathology/Braak NFT and
                                      Celltype.l1/l2/l3 labels. Barcodes are "LibraryNN_<10x barcode>", NN = position
                                      of the library in L1..L8, K1..K8, J1..J8 (checked: first barcodes match).
  GSE295998_family.soft.gz            per-library subject, sex, age, disease state, cohort, pathology, Braak NFT.
Counts: CellBender "filtered" matrices; values are rounded to integers (np.rint) in case CellBender wrote floats
(in the files checked they were already integer-valued). Cells kept = barcodes present in the authors' cluster file
(i.e. the authors' QC + doublet removal), then the uniform QC (>= 200 genes, < 10% mito UMIs) on the CellBender counts.
author_cell_type = Celltype.l3 (finest; Celltype.l2 where l3 is NA: 93 ciliated-glia cells); Celltype.l1/l2 kept as extra columns.
donor = Subject (Subject01 contributes two libraries: L1 and K1). disease: PD / control (no dementia info -> PD).
braak_lb: NaN (the GEO "Braak stage" is Braak & Braak NFT -> braak_nft); Lewy pathology type in lb_pathology.
"""
import gzip
import re
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from pdrep.harmonise import RAW

D = RAW / "glock2026"
LIBS = [f"{p}{i}" for p in "LKJ" for i in range(1, 9)]  # Library01..24 order

AUTHOR_TO_CLASS = {}
for _n in range(1, 12):
    AUTHOR_TO_CLASS.update({f"Oligo.{_n}": "Oligodendrocytes", f"Ex.{_n}": "Excitatory neurons",
                            f"In.{_n}": "Inhibitory neurons", f"MG/PVM.{_n}": "Microglia",
                            f"Astro.{_n}": "Astrocytes", f"OPC.{_n}": "OPCs", f"CG.{_n}": "Ependymal"})
AUTHOR_TO_CLASS.update({"COP": "OPCs", "Capillary": "Endothelial", "Venous": "Endothelial", "Arterial": "Endothelial",
                        "Pericyte": "Pericytes", "SMC": "Pericytes", "PV_Fibroblast": "Fibroblasts", "Tcell": "T cells",
                        "Bcell": "Immune (other)", "Precursor": "Neurons (other)",
                        # "Mix": Celltype.l2 = Excitatory, subtype unclear; "Ciliated_glia": l2 fallback for 93 CG cells
                        # without an l3 label; MG/PVM.8 and .10 are perivascular-macrophage-like (MRC1/F13A1/CD163 high,
                        # P2RY12 low) -> Immune (other)
                        "Mix": "Neurons (other)", "Ciliated_glia": "Ependymal",
                        "MG/PVM.8": "Immune (other)", "MG/PVM.10": "Immune (other)"})


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
            elif cur is not None and line.startswith("!Sample_description"):
                cur["library"] = line.split("=", 1)[1].split(":")[-1].strip()
            elif cur is not None and line.startswith("!Sample_characteristics_ch1"):
                k, _, v = line.split("=", 1)[1].strip().partition(":")
                cur[k.strip().lower()] = v.strip()
    return pd.DataFrame(rows).set_index("library")


def read_cb(path):
    with h5py.File(path, "r") as f:
        m = f["matrix"]
        X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=tuple(m["shape"][:])).T.tocsr()
        bc = m["barcodes"][:].astype(str)
        ft = m["features"]
        ftype = ft["feature_type"][:].astype(str)
        var = pd.DataFrame({"gene_ids": ft["id"][:].astype(str)}, index=ft["name"][:].astype(str))
    keep = ftype == "Gene Expression"
    nonint = float(np.mean(X.data != np.rint(X.data))) if X.nnz else 0.0
    X.data = np.rint(X.data)
    X.eliminate_zeros()
    a = ad.AnnData(X[:, keep], obs=pd.DataFrame(index=bc), var=var[keep])
    a.var_names_make_unique()
    return a, nonint


def qc(a, min_genes=200, max_mito=10.0):
    X = a.X.tocsr()
    ng = np.diff(X.indptr)
    tot = np.asarray(X.sum(1)).ravel()
    mt = np.asarray(X[:, a.var_names.str.startswith("MT-")].sum(1)).ravel()
    keep = (ng >= min_genes) & (100 * mt / np.maximum(tot, 1) < max_mito)
    return a[keep].copy()


def load():
    meta = soft_samples(D / "GSE295998_family.soft.gz")
    u = pd.read_csv(D / "GSE295998_UmapCoordinatesAndClusterAssignments.txt.gz", sep="\t",
                    usecols=["Barcode", "Library", "Subject", "Condition", "Pathology", "Braak_braak_nft_stage",
                             "Celltype.l1", "Celltype.l2", "Celltype.l3"]).set_index("Barcode")
    parts = []
    for f in sorted(D.glob("GSM*_cellbender_filtered.h5")):
        gsm, lib = re.match(r"(GSM\d+)_(\w+?)_cellbender_filtered\.h5", f.name).groups()
        a, nonint = read_cb(f)
        n0 = a.n_obs
        a.obs_names = f"Library{LIBS.index(lib) + 1:02d}_" + a.obs_names
        ul = u[u.Library == lib]
        a = a[a.obs_names.isin(ul.index)].copy()
        n1 = a.n_obs
        a = qc(a)
        print(f"{lib}: cellbender {n0}, in author table {n1} (author table {len(ul)}), after QC {a.n_obs}; "
              f"non-integer fraction {nonint:.4f}", flush=True)
        m = meta.loc[lib]
        ua = ul.loc[a.obs_names]
        a.obs["sample"] = lib
        a.obs["gsm"] = gsm
        a.obs["donor"] = ua["Subject"].values
        a.obs["region"] = "amygdala"
        a.obs["disease"] = "PD" if m["disease state"] == "PD" else "control"
        a.obs["disease_raw"] = m["disease state"]
        a.obs["lb_pathology"] = m["pathology"]
        a.obs["sex"] = {"Female": "F", "Male": "M"}[m["sex"]]
        a.obs["age"] = float(m["age"])
        a.obs["pmi"] = np.nan
        a.obs["braak_lb"] = np.nan
        a.obs["braak_nft"] = float(re.search(r"(\d)", m["braak and_braak_nft_stage"]).group(1))
        a.obs["cohort"] = m["cohort"]
        a.obs["author_cell_type"] = ua["Celltype.l3"].fillna(ua["Celltype.l2"]).values  # l3, else l2
        a.obs["author_l1"] = ua["Celltype.l1"].values
        a.obs["author_l2"] = ua["Celltype.l2"].values
        assert (ua["Condition"] == m["disease state"]).all()
        parts.append(a)
    return ad.concat(parts, join="inner", merge="same")

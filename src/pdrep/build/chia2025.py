"""chia2025 — Chia et al. 2025 Genome Med (GSE303823), 10x snRNA-seq of prefrontal cortex (middle frontal gyrus):
4 NCI controls, 4 AD, 6 DLB, 6 PDD (one library per donor).

Sources (raw/chia2025/):
  GSE303823_rna_annotated.rds: Seurat v4 object. Converted on SLURM by scripts/convert/chia2025.R (env mast_rev, no
      SeuratObject; slots read via attr()) into converted/: counts.mtx.gz (RNA assay @counts: raw integer UMI,
      genes x cells), genes.tsv, barcodes.tsv, meta.tsv (meta.data: sample, group, sex, cell.type, ...).
      SCT and data slots were not used.
  raw/supplements/chia2025/13073_2025_1519_MOESM1_ESM.xlsx: Table S1 (per-sample B&B stage and cell counts per type,
      under IDs NCI_1-4/AD_1-4/DLB_1-7/PDD_1-6) and Table S3 (neuropathology). The GEO sample IDs (CTRL_4-7, ADD_4-7,
      DLB_2/4/5/6/7/8, PDD_1/2/4/5/7/9) are matched to the table IDs in order within each group; this was verified by
      identical per-cell-type counts (e.g. ADD_4 = AD_1: 1956 astrocytes, 3766 ExN). Age is only given as group mean.
"""
import numpy as np
import pandas as pd
import scipy.io
import anndata as ad
from pdrep.harmonise import RAW

D = RAW / "chia2025" / "converted"
AUTHOR_TO_CLASS = {"Oligodendrocyte": "Oligodendrocytes", "Inhibitory neuron": "Inhibitory neurons",
                   "Excitatory neuron": "Excitatory neurons", "Astrocyte": "Astrocytes", "OPC": "OPCs",
                   "Microglia": "Microglia", "Ependymal cell": "Ependymal", "Endothelial cell": "Endothelial"}
DISEASE = {"CTRL": "control", "ADD": "AD", "DLB": "DLB", "PDD": "PDD"}
SUPPL_ID = {"CTRL_4": "NCI_1", "CTRL_5": "NCI_2", "CTRL_6": "NCI_3", "CTRL_7": "NCI_4",
            "ADD_4": "AD_1", "ADD_5": "AD_2", "ADD_6": "AD_3", "ADD_7": "AD_4",
            "DLB_2": "DLB_1", "DLB_4": "DLB_2", "DLB_5": "DLB_3", "DLB_6": "DLB_4", "DLB_7": "DLB_5", "DLB_8": "DLB_7",
            "PDD_1": "PDD_1", "PDD_2": "PDD_2", "PDD_4": "PDD_3", "PDD_5": "PDD_4", "PDD_7": "PDD_5", "PDD_9": "PDD_6"}
# Table S1 "B&B stage" (Braak & Braak NFT stage; "-" = not given)
BRAAK_NFT = {"NCI_1": "III", "NCI_2": "II", "NCI_3": "III", "NCI_4": "I", "AD_1": "VI", "AD_2": "VI", "AD_3": "VI",
             "AD_4": "VI", "DLB_4": "III", "DLB_7": "V", "PDD_6": "II"}


def load():
    X = scipy.io.mmread(D / "counts.mtx.gz").T.tocsr()
    genes = pd.read_csv(D / "genes.tsv", header=None)[0].astype(str).values
    cells = pd.read_csv(D / "barcodes.tsv", header=None)[0].astype(str).values
    m = pd.read_csv(D / "meta.tsv", sep="\t", index_col=0).loc[cells]
    sid = m["sample"].map(SUPPL_ID)
    assert sid.notna().all()
    obs = pd.DataFrame({
        "donor": m["sample"].values, "sample": m["sample"].values, "suppl_id": sid.values,
        "region": "PFC", "disease": m.group.map(DISEASE).values, "disease_raw": m.group.values,
        "sex": m.sex.values, "age": np.nan, "pmi": np.nan, "braak_lb": np.nan,
        "braak_nft": sid.map(BRAAK_NFT).values,
        "author_cell_type": m["cell.type"].values, "seurat_cluster_author": m.seurat_clusters.astype(str).values,
    }, index=cells)
    assert obs.disease.notna().all()
    return ad.AnnData(X.astype(np.float32), obs=obs, var=pd.DataFrame(index=genes))

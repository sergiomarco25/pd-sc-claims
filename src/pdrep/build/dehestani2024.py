"""dehestani2024 — Dehestani et al. 2024 Mol Brain (GSE272760), 10x 3' v3.1 snRNA-seq of PFC and ACC from
2 LRRK2-PD, 2 GBA1-PD and 2 control donors (12 libraries = 6 donors x 2 regions).

Sources (raw/dehestani2024/):
  GSE272760_all_new_SCT_harmony_res0.1_v1.rds.gz: Seurat v4 object (gzip-compressed twice). Converted on SLURM by
      scripts/convert/dehestani2024.R (env mast_rev, no SeuratObject; slots read via attr()) into converted/:
      counts.mtx.gz (RNA assay @counts: raw integer UMI, genes x cells), genes.tsv, barcodes.tsv, meta.tsv (meta.data).
      SCT and data slots were not used.
  raw/supplements/dehestani2024/13041_2024_1128_MOESM1_ESM.xlsx (Table S1): per-sample mutation, region, Braak LB
      stage, sex, age, brain bank code (= donor), postmortem delay, pH.
"""
import numpy as np
import pandas as pd
import scipy.io
import anndata as ad
from pdrep.harmonise import RAW

D = RAW / "dehestani2024" / "converted"
S1 = RAW / "supplements" / "dehestani2024" / "13041_2024_1128_MOESM1_ESM.xlsx"
AUTHOR_TO_CLASS = {"ODC": "Oligodendrocytes", "Astro": "Astrocytes", "MG": "Microglia", "OPC": "OPCs",
                   **{f"ExN{i}": "Excitatory neurons" for i in range(1, 7)},
                   **{f"InN{i}": "Inhibitory neurons" for i in range(1, 3)},
                   "Vas": "Endothelial"}


def _pmi(v):
    v = str(v)
    if ":" in v:
        h, m = v.split(":")[:2]
        return int(h) + int(m) / 60
    try:
        return float(v)
    except ValueError:
        return np.nan


def load():
    X = scipy.io.mmread(D / "counts.mtx.gz").T.tocsr()
    genes = pd.read_csv(D / "genes.tsv", header=None)[0].astype(str).values
    cells = pd.read_csv(D / "barcodes.tsv", header=None)[0].astype(str).values
    m = pd.read_csv(D / "meta.tsv", sep="\t", index_col=0).loc[cells]
    s = pd.read_excel(S1, header=None).iloc[3:15, [0, 1, 2, 3, 4, 5, 7, 8, 9]]
    s.columns = ["sample", "mutation", "region", "braak", "sex", "age", "bank_code", "pmi", "ph"]
    s["sample"] = s["sample"].str.replace("sample", "Sample")
    s = s.set_index("sample")
    r = s.loc[m.SamplID.values]
    assert (r.region.values == m.Region.values).all()
    geno = r.mutation.astype(str).str.split(" ").str[0]
    assert (geno.str.replace("GBA1", "GBA").values == m.Mutation.values).all()
    obs = pd.DataFrame({
        "donor": r.bank_code.values, "sample": m.SamplID.values,
        "region": m.Region.map({"PFC": "PFC", "ACC": "cingulate"}).values,
        "region_author": m.Region.values,
        "disease": np.where(geno.values == "HC", "control", "PD"),
        "disease_raw": np.where(geno.values == "HC", "HC", "PD " + r.mutation.astype(str).values),
        "genotype": np.where(geno.values == "HC", "none", geno.values),
        "mutation": r.mutation.astype(str).values,
        "sex": r.sex.values, "age": r.age.astype(float).values, "pmi": [_pmi(v) for v in r.pmi],
        "braak_lb": np.nan, "braak_lb_group": r.braak.astype(str).values, "brain_ph": r.ph.astype(float).values,
        "author_cell_type": m.CellType.values, "author_cell_type_broad": m.BroadCellType.values,
    }, index=cells)
    return ad.AnnData(X.astype(np.float32), obs=obs, var=pd.DataFrame(index=genes))

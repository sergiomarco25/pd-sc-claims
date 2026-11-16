"""smajic2022 — Smajić et al. 2022 Brain, GSE157783: ventral midbrain snRNA-seq, 5 IPD / 6 control donors.

Adapted from revision/scripts/06a_external_build_and_score.py::load_smajic2022_midbrain (logic copied, not imported).
Counts: GSE157783_IPDCO_hg_midbrain_UMI.tsv (dense genes x cells UMI table, 41,435 post-QC nuclei), gene rows from
_genes.tsv (Ensembl IDs, mapped to symbols with the GRCh38 10x features file of GSE161045), author labels and donor
from _cell.tsv (cell_ontology, patient). The cell table codes IPD donors as PD1..PD5; the nuclei per patient match
Supplementary Table 1 "Profiled nuclei" for IPD1..IPD5 exactly, so PDk = IPDk (GEO title).
Donor metadata: GEO soft (sex, age, PMI) + Supplementary Table 1 (Supplementary data 1 PDF, p.7: LB Braak stage,
onset, duration, non-motor features, SN neuron loss), transcribed below.
"""
import gzip
import tarfile
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

from pdrep import config as cfg  # noqa: E402

DATA = cfg.DATA_ROOT
GEO = DATA / "external_validation"
EXT = GEO / "_extracted"
ACC = "GSE157783"

AUTHOR_TO_CLASS = {
    "Astrocytes": "Astrocytes", "CADPS2+ neurons": "Neurons (other)", "DaNs": "DA neurons",
    "Endothelial cells": "Endothelial", "Ependymal": "Ependymal", "Excitatory": "Excitatory neurons",
    "GABA": "Inhibitory neurons", "Inhibitory": "Inhibitory neurons", "Microglia": "Microglia",
    "Oligodendrocytes": "Oligodendrocytes", "OPCs": "OPCs", "Pericytes": "Pericytes",
}

# Supplementary Table 1 (single-cell cohort): onset, duration, LB Braak, non-motor features, SN neuron loss.
SUPP_T1 = pd.DataFrame([
    ("IPD1", 71, 13, 5, "Depression, visual hallucinations", "severe"),
    ("IPD2", 57, 10, 6, "Hallucinations, paranoia", "severe"),
    ("IPD3", 71, 7, 6, "Stroke, cognitive impairment, dementia", "severe"),
    ("IPD4", 69, 12, 6, "Cognitive impairment, hallucinations, depression", "severe"),
    ("IPD5", 71, 8, 5, "Blindness, memory impairment", "severe"),
    ("C1", np.nan, np.nan, np.nan, "", "mild"),
    ("C2", np.nan, np.nan, np.nan, "", "none"),
    ("C3", np.nan, np.nan, np.nan, "", "none"),
    ("C4", np.nan, np.nan, np.nan, "", "none"),
    ("C5", np.nan, np.nan, np.nan, "", "none"),
    ("C6", np.nan, np.nan, np.nan, "", "none"),
], columns=["donor", "age_onset", "disease_duration", "braak_lb", "nonmotor_features", "sn_neuron_loss"]).set_index("donor")


def soft(acc):
    rows, cur = [], None
    for l in gzip.open(GEO / acc / f"{acc}_family.soft.gz", "rt"):
        if l.startswith("^SAMPLE"):
            cur = {"gsm": l.split("=")[1].strip()}; rows.append(cur)
        elif cur is not None and l.startswith("!Sample_title"):
            cur["title"] = l.split("=", 1)[1].strip()
        elif cur is not None and l.startswith("!Sample_characteristics_ch1"):
            kv = l.split("=", 1)[1].strip().split(":", 1)
            if len(kv) == 2:
                cur[kv[0].strip()] = kv[1].strip()
    return pd.DataFrame(rows)


def _tsv(name):
    t = tarfile.open(GEO / ACC / f"{ACC}_IPDCO_hg_midbrain_{name}.tar.gz")
    return t.extractfile(t.getmembers()[0])


def load():
    cells = pd.read_csv(_tsv("cell"), sep="\t")
    ens = pd.read_csv(_tsv("genes"), sep="\t")["gene"].astype(str).values      # Ensembl IDs, row order of UMI table
    fx = next((EXT / "GSE161045").glob("*features.tsv.gz"))                  # Ensembl -> symbol (GRCh38, 10x)
    f2s = pd.read_csv(fx, sep="\t", header=None); f2s[0] = f2s[0].str.split(".").str[0]
    f2s = f2s.drop_duplicates(0).set_index(0)[1]
    chunks, cols = [], None
    for ch in pd.read_csv(_tsv("UMI"), sep="\t", chunksize=2000):           # header = cells, rows = genes
        cols = ch.columns
        chunks.append(sp.csr_matrix(ch.values.astype(np.float32)))
    M = sp.vstack(chunks).T.tocsr()                                          # cells x genes
    assert M.shape[1] == len(ens), (M.shape, len(ens))
    a = ad.AnnData(M, obs=pd.DataFrame(index=cols.astype(str)))
    a.var["gene_ids"] = ens
    a.var_names = [f2s.get(e, e) for e in ens]
    a.var_names_make_unique()

    cells = cells.set_index("barcode").reindex(a.obs_names)
    assert cells["patient"].notna().all()
    s = soft(ACC).set_index("title")
    donor = cells["patient"].astype(str).str.replace(r"^PD", "IPD", regex=True)
    a.obs["donor"] = donor.values
    a.obs["sample"] = a.obs["donor"]                       # one library per donor (barcode suffix _1.._11)
    a.obs["library_index"] = a.obs_names.str.split("_").str[-1].values
    a.obs["gsm"] = donor.map(s["gsm"]).values
    a.obs["region"] = "midbrain"
    a.obs["disease_raw"] = donor.map(s["disease state"]).values
    t1 = SUPP_T1.reindex(donor.values)
    dementia = t1["nonmotor_features"].str.contains("dementia", case=False).values
    a.obs["disease"] = np.where(a.obs.disease_raw.str.startswith("Idiopathic"), np.where(dementia, "PDD", "PD"),
                                "control")
    a.obs["sex"] = donor.map(s["Sex"]).values
    a.obs["age"] = donor.map(s["age at death"]).astype(float).values
    a.obs["pmi"] = donor.map(s["pmi (h)"]).astype(float).values
    for c in SUPP_T1.columns:
        a.obs[c] = t1[c].values
    a.obs["author_cell_type"] = cells["cell_ontology"].astype(str).values
    assert a.obs.donor.nunique() == 11 and a.obs.gsm.notna().all()
    return a

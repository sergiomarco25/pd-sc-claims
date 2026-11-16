"""kamath2022 — Kamath et al. 2022 Nat Neurosci, GSE178265 (human part): SNpc snRNA-seq of 8 control + 10 PD/LBD
donors (NR4A2-sorted DA-enriched and unsorted DAPI libraries) and dorsal striatum (caudate) of 4 control donors.

Adapted from revision/scripts/06a_external_build_and_score.py::load_kamath2022_sn (logic copied, not imported).
Changes: keeps ALL human libraries (the revision code kept only FACS-negative SN libraries), adds `enrichment`
(sorted = NR4A2-positive "Pos*" or NeuN-sorted; unsorted = DAPI / NR4A2-negative), full donor metadata.
Counts: GSE178265_Homo_matrix.mtx.gz (genes x cells, 434,340 post-QC nuclei; hg19 gene symbols only, no Ensembl IDs),
barcodes prefixed by library name = GEO sample title. Streamed with pyarrow (5 GB gz).
Donor metadata: GEO soft (tissue, sex, age, PMI, disease, status, cause of death, FACS classification).
Library pPDCN5730NeuN22119 is annotated "Substantia nigra pars compacta" on GEO but is a caudate (CN) NeuN-sorted
library of a caudate-only donor; region is taken from the title (CN -> caudate).
No author cell labels on disk (GEO deposit has none; Broad SCP metadata not downloaded) -> AUTHOR_TO_CLASS = {}.
"""
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

from pdrep import config as cfg  # noqa: E402

DATA = cfg.DATA_ROOT
GEO = DATA / "external_validation"
ACC = "GSE178265"
AUTHOR_TO_CLASS = {}


def soft(acc):
    rows, cur = [], None
    for l in gzip.open(GEO / acc / f"{acc}_family.soft.gz", "rt"):
        if l.startswith("^SAMPLE"):
            cur = {"gsm": l.split("=")[1].strip()}; rows.append(cur)
        elif cur is not None and l.startswith("!Sample_title"):
            cur["title"] = l.split("=", 1)[1].strip()
        elif cur is not None and l.startswith("!Sample_organism_ch1"):
            cur["organism"] = l.split("=", 1)[1].strip()
        elif cur is not None and l.startswith("!Sample_characteristics_ch1"):
            kv = l.split("=", 1)[1].strip().split(":", 1)
            if len(kv) == 2:
                cur[kv[0].strip()] = kv[1].strip()
    return pd.DataFrame(rows)


def sample_meta():
    s = soft(ACC)
    s = s[s.organism == "Homo sapiens"].copy()
    s["donor"] = s.title.str.extract(r"(?:CN|SNxi)(\d{4})")[0]
    s["region"] = np.where(s.title.str.startswith("pPDCN"), "caudate", "SN")
    s["enrichment"] = np.where(s["facs classification"].isin(["Positive", "NeuN"]), "sorted", "unsorted")
    s["disease_h"] = s["status"].map({"Ctrl": "control", "PD": "PD", "LBD": "LBD"})
    assert s.donor.notna().all() and s.disease_h.notna().all()
    return s.set_index("title")


def load():
    meta = sample_meta()
    bcd = pd.read_csv(GEO / ACC / f"{ACC}_Homo_bcd.tsv.gz", header=None)[0].astype(str)
    feats = pd.read_csv(GEO / ACC / f"{ACC}_Homo_features.tsv.gz", header=None, sep="\t")[0].astype(str).values
    samp = bcd.str.rsplit("_", n=1).str[0]
    assert samp.isin(meta.index).all(), set(samp) - set(meta.index)
    import pyarrow.csv as pacsv
    opts = pacsv.ReadOptions(skip_rows=2, column_names=["g", "c", "v"], block_size=1 << 28)
    reader = pacsv.open_csv(GEO / ACC / f"{ACC}_Homo_matrix.mtx.gz", read_options=opts,
                            parse_options=pacsv.ParseOptions(delimiter=" "))
    rows, cols, vals = [], [], []
    for batch in reader:
        rows.append((batch.column("c").to_numpy() - 1).astype(np.int32))
        cols.append((batch.column("g").to_numpy() - 1).astype(np.int32))
        vals.append(batch.column("v").to_numpy().astype(np.float32))
    r = np.concatenate(rows); del rows
    c = np.concatenate(cols); del cols
    v = np.concatenate(vals); del vals
    M = sp.csr_matrix((v, (r, c)), shape=(len(bcd), len(feats)))
    del r, c, v
    a = ad.AnnData(M, obs=pd.DataFrame(index=bcd.values))
    a.var_names = feats
    a.var_names_make_unique()
    m = meta.loc[samp.values]
    a.obs["donor"] = m["donor"].values
    a.obs["sample"] = samp.values
    a.obs["gsm"] = m["gsm"].values
    a.obs["region"] = m["region"].values
    a.obs["disease"] = m["disease_h"].values
    a.obs["disease_raw"] = (m["disease"] + " / " + m["status"]).values
    a.obs["sex"] = m["Sex"].map({"Male": "M", "Female": "F"}).values
    a.obs["age"] = m["age"].astype(float).values
    a.obs["pmi"] = m["pmi"].astype(float).values
    a.obs["braak_lb"] = np.nan
    a.obs["cause_of_death"] = m["Cause of death"].values
    a.obs["facs_classification"] = m["facs classification"].values
    a.obs["enrichment"] = m["enrichment"].values
    a.obs["author_cell_type"] = np.nan
    return a

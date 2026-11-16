"""xu2023 — Xu et al. 2023 Nat Aging, GSE161045: putamen snRNA-seq, 4 control / 4 AD / 4 PD.

Adapted from revision/scripts/06a_external_build_and_score.py::load_xu2023_putamen (logic copied, not imported).
Counts: per-sample 10x matrices from GSE161045_RAW.tar (already extracted in external_validation/_extracted/GSE161045);
all deposited barcodes are kept (38,929 = the paper's pre-QC count; the paper's 30,908 post-QC barcodes are not deposited).
Donor metadata: Supplementary Table 1 (MOESM3 xlsx, sheet Table_S1_Patient_information).
No author cell labels deposited -> AUTHOR_TO_CLASS = {} (marker clustering in finalise).
"""
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
import anndata as ad

from pdrep import config as cfg  # noqa: E402

DATA = cfg.DATA_ROOT
GEO = DATA / "external_validation"
EXT = GEO / "_extracted"
SUPP = DATA / "pd_sc_replication" / "raw" / "supplements" / "xu2023" / "43587_2023_363_MOESM3_ESM.xlsx"
ACC = "GSE161045"
AUTHOR_TO_CLASS = {}


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


def read_10x_prefixed(folder, prefix):
    """10x triplet stored as <prefix>_barcodes / features|genes / matrix files (read without symlinking)."""
    m = scipy.io.mmread(folder / f"{prefix}_matrix.mtx.gz").T.tocsr()
    bc = pd.read_csv(folder / f"{prefix}_barcodes.tsv.gz", header=None)[0].astype(str).values
    ff = folder / f"{prefix}_features.tsv.gz"
    feat = pd.read_csv(ff if ff.exists() else folder / f"{prefix}_genes.tsv.gz", header=None, sep="\t")
    var = pd.DataFrame({"gene_ids": feat[0].astype(str).values}, index=feat[1].astype(str).values)
    a = ad.AnnData(sp.csr_matrix(m, dtype=np.float32), obs=pd.DataFrame(index=bc), var=var)
    a.var_names_make_unique()
    return a


def donor_table():
    t = pd.read_excel(SUPP, sheet_name="Table_S1_Patient_information", header=2)
    t.columns = [str(c).strip() for c in t.columns]
    t = t[t["Case"].notna() & t["Case"].astype(str).str.contains("-")].copy()
    t["Case"] = t["Case"].astype(str).str.strip()
    return t.set_index("Case")


def load():
    raw = EXT / ACC
    s = soft(ACC)
    meta = donor_table()
    ads = []
    for _, r in s.iterrows():
        prefix = next(raw.glob(r.gsm + "_*_matrix.mtx.gz")).name.replace("_matrix.mtx.gz", "")
        a = read_10x_prefixed(raw, prefix)
        short = r.title.replace("Putamen_", "")                 # CTRL1 / AD1 / PD1
        donor = short.rstrip("0123456789") + "-" + short[len(short.rstrip("0123456789")):]  # CTRL-1 (Table S1)
        m = meta.loc[donor]
        a.obs["donor"] = donor
        a.obs["sample"] = r.title
        a.obs["gsm"] = r.gsm
        a.obs["region"] = "putamen"
        a.obs["disease_raw"] = r.diagnosis
        a.obs["disease"] = {"control": "control", "Alzheimer disease (AD)": "AD",
                            "Parkinson disease (PD)": "PD"}[r.diagnosis]
        a.obs["sex"] = str(m["Sex"]).strip()
        a.obs["age"] = float(m["Age at death (years)"])
        a.obs["pmi"] = float(m["PMI:  (hours)"])
        a.obs["braak_lb"] = float(m["Braak Lewy body stage"])
        a.obs["braak_nft"] = str(m["Braak NFT stage"]).strip()
        a.obs["braak_abeta"] = str(m["Braak amyloid-beta stage"]).strip()
        a.obs["cdr"] = float(m["Clinical Dementia Rating"])
        a.obs["cause_of_death"] = str(m["Cause of Death"]).strip()
        a.obs["author_cell_type"] = np.nan
        a.obs_names = r.gsm + "_" + a.obs_names
        ads.append(a)
    out = ad.concat(ads, merge="same")
    assert out.obs.donor.nunique() == 12
    return out

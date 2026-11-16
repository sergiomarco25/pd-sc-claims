"""zhu2024 — Zhu et al. 2024 Sci Transl Med, GSE202210: prefrontal cortex (BA9) snRNA-seq, 6 PD / 6 control.

Adapted from revision/scripts/06a_external_build_and_score.py::load_zhu2024_pfc (logic copied, not imported).
Counts: per-sample 10x matrices from GSE202210_RAW.tar (already extracted in external_validation/_extracted/GSE202210).
Donor metadata: GEO soft gives only diagnosis; Table S1 (sex/age/PMI/LB score) is not on disk
(raw/supplements/zhu2024 is empty) -> sex/age/pmi/braak_lb NaN.
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
ACC = "GSE202210"
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


def load():
    raw = EXT / ACC
    s = soft(ACC)
    ads = []
    for _, r in s.iterrows():
        prefix = next(raw.glob(r.gsm + "_*_matrix.mtx.gz")).name.replace("_matrix.mtx.gz", "")
        a = read_10x_prefixed(raw, prefix)
        donor = r.title.replace("_HHT", "")                  # e.g. hsDG142PD (one library per donor)
        a.obs["donor"] = donor
        a.obs["sample"] = r.gsm
        a.obs["region"] = "PFC"
        a.obs["disease_raw"] = r.diagnosis
        a.obs["disease"] = "PD" if "Parkinson" in r.diagnosis else ("control" if r.diagnosis == "Control" else "other")
        for c in ("sex", "age", "pmi", "braak_lb", "author_cell_type"):
            a.obs[c] = np.nan
        a.obs_names = r.gsm + "_" + a.obs_names
        ads.append(a)
    out = ad.concat(ads, merge="same")
    assert out.obs.donor.nunique() == 12
    return out

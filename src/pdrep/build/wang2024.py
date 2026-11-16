"""wang2024 — Wang et al. 2024 Sci Adv, GSE184950: substantia nigra snRNA-seq, 23 PD (6 PD + 17 PDD) / 9 control donors.

Adapted from revision/scripts/06a_external_build_and_score.py::load_wang2024_sn (logic copied, not imported).
Counts: per-library Cell Ranger filtered_feature_bc_matrix from GSE184950_RAW.tar (already extracted in
external_validation/_extracted/GSE184950). 34 GEO libraries = 32 donors (HBNX: p7+p9, HBHF: p4+p5).
Donor metadata: GEO soft characteristics (donor id, age, sex, PMI, Braak stage, diagnosis), cross-checked with
Table S1 (raw/supplements/wang2024, sheet "S1. Sample cohort").
QC: the deposit holds all 457,453 pre-QC nuclei (paper: "457,453 ... before QC"). The paper's first QC step is applied
here: keep nuclei with 200 <= detected genes <= 2500 (paper: 355,157 retained). The later removals (doublets /
unstable clusters, down to 315,867) are not reproducible from the deposit. The "Braak stage" (given for PD only, N/A for all
controls) is taken as the Lewy body (PD) Braak stage.
No author cell labels deposited (GSE184950_add2.xlsx is only the GEO metadata template) -> AUTHOR_TO_CLASS = {}.
"""
import gzip
from pathlib import Path
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc

from pdrep import config as cfg  # noqa: E402

DATA = cfg.DATA_ROOT
GEO = DATA / "external_validation"
EXT = GEO / "_extracted"
ACC = "GSE184950"
AUTHOR_TO_CLASS = {}
MIN_GENES, MAX_GENES = 200, 2500             # paper Methods, "snRNA-seq data preprocessing"
ROMAN = {"Stage I": 1, "Stage II": 2, "Stage III": 3, "Stage IV": 4, "Stage V": 5, "Stage VI": 6}


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


def load():
    raw = EXT / ACC
    s = soft(ACC)
    ads = []
    for _, r in s.iterrows():
        dd = next(p for p in raw.glob(r.gsm + "_*") if p.is_dir())
        mdir = next(p for p in dd.rglob("matrix.mtx.gz")).parent
        a = sc.read_10x_mtx(mdir, var_names="gene_symbols", make_unique=True)
        ds = r["disease state"]
        a.obs["donor"] = r["brain bank donor id"]
        a.obs["sample"] = r.title
        a.obs["gsm"] = r.gsm
        a.obs["region"] = "SN"
        a.obs["disease_raw"] = ds
        a.obs["disease"] = {"Unaffected Control": "control", "Parkinson's Disease": "PD",
                            "Parkinson's Disease Dementia": "PDD"}[ds]
        a.obs["sex"] = {"Male": "M", "Female": "F"}[r.gender]
        a.obs["age"] = float(r.age)
        a.obs["pmi"] = float(r["postmortem interval hours"])
        a.obs["braak_lb"] = ROMAN.get(r["braak stage"], np.nan)
        a.obs["braak_stage_raw"] = r["braak stage"]
        a.obs["race"] = r.race
        a.obs["author_cell_type"] = np.nan
        a.obs_names = r.gsm + "_" + a.obs_names
        ng = np.asarray((a.X > 0).sum(axis=1)).ravel()
        a = a[(ng >= MIN_GENES) & (ng <= MAX_GENES)].copy()
        ads.append(a)
    out = ad.concat(ads, merge="same")
    assert out.obs.donor.nunique() == 32
    return out

"""ma2025 — Ma et al. 2025 Nat Commun (GSE253462), snRNA-seq of PD substantia nigra and cingulate cortex.

Sources (raw/ma2025/):
  GSE253462_RAW.tar -> extracted/GSM*_<seqid>_{nigra,cingulate}.txt.gz: one count table per sample, genes x cells,
      space-separated R write.table (header = cell barcodes, first column = gene symbol). Integer counts after
      ambient removal (CellBender for SN, DecontX for cingulate; authors' QC-filtered nuclei). Cingulate tables have
      50,360 genes, nigra tables 36,601 (different references); concatenated with an outer join, absent genes = 0.
  GSE253462_family.soft.gz: per-sample donor (TB number), tissue, sex, exact age, condition.
  raw/supplements/ma2025/41467_2025_62478_MOESM3_ESM.xlsx sheet A (Supplementary Data 1): neuropath diagnosis
      (PD vs PD with dementia), Braak NFT, CERAD per donor x sequencing ID.
No per-cell author labels are public (only per-sample counts per type, MOESM3 sheet C) -> marker clustering.
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from pdrep.harmonise import RAW

AUTHOR_TO_CLASS = {}
D = RAW / "ma2025"
SUPP = RAW / "supplements" / "ma2025" / "41467_2025_62478_MOESM3_ESM.xlsx"


def _soft():
    rows, cur = [], None
    import gzip
    with gzip.open(D / "GSE253462_family.soft.gz", "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith("^SAMPLE"):
                cur = {"gsm": line.split(" = ")[1]}
                rows.append(cur)
            elif cur is not None and line.startswith("!Sample_characteristics_ch1"):
                k, v = line.split(" = ", 1)[1].split(": ", 1)
                cur[k.strip().lower()] = v.strip()
            elif cur is not None and line.startswith("!Sample_supplementary_file_1"):
                cur["file"] = line.split("/")[-1]
            elif cur is not None and line.startswith("!Sample_title"):
                cur["title"] = line.split(" = ")[1].split(",")[0]
    s = pd.DataFrame(rows)
    s["seq_id"] = s.file.str.extract(r"GSM\d+_([A-Z]+\d+)_")[0]
    return s


def _read_table(path, chunksize=4000):
    """genes x cells space-separated table -> (cells x genes csr, genes, barcodes)."""
    with pd.read_csv(path, sep=" ", index_col=0, chunksize=chunksize) as it:
        mats, genes = [], []
        for ch in it:
            v = ch.to_numpy()
            assert np.array_equal(v, np.round(v)), f"{path.name}: non-integer counts"
            mats.append(sp.csr_matrix(v.astype(np.float32)))
            genes.extend(ch.index.astype(str))
            cells = ch.columns.astype(str)
    return sp.vstack(mats).T.tocsr(), genes, list(cells)


def load():
    s = _soft()
    meta = pd.read_excel(SUPP, sheet_name="A", header=0).iloc[:, :13]
    meta.columns = ["tcr", "snrna", "st", "tb", "seq_id", "brain_region", "condition", "neuropath_dx", "age_bin",
                    "sex", "braak_nft", "path_braak_asyn", "cerad"]
    meta["seq_num"] = meta.seq_id.str.replace(r"^(OF|OJ)0*", "", regex=True)
    parts = []
    for _, r in s.iterrows():
        X, genes, cells = _read_table(D / "extracted" / r.file)
        reg = "SN" if "nigra" in r.file else "cingulate"
        prefix = "OJ" if reg == "SN" else "OF"
        num = re.sub(r"^(OJ|T)", "", r.seq_id)
        m = meta[meta.seq_id.str.startswith(prefix) & (meta.seq_num == num)] if reg == "SN" else \
            meta[(meta.tb.str.replace("-", "") == r.seq_id) & (meta.brain_region == "Cingulate Cortex")]
        assert len(m) == 1, (r.file, len(m))
        m = m.iloc[0]
        a = ad.AnnData(X, var=pd.DataFrame(index=genes))
        a.obs_names = [f"{r.seq_id}_{c}" for c in cells]
        age = r.age.rstrip("+")
        dx = str(m.neuropath_dx)
        a.obs = a.obs.assign(
            donor=m.tb, sample=m.seq_id, region=reg,
            disease=("PDD" if "with dementia" in dx else "PD") if r.condition == "PD" else "control",
            disease_raw=f"{r.condition}; {dx}", sex=r.sex, age=float(age), age_censored_89=r.age.endswith("+"),
            pmi=np.nan, braak_lb=np.nan, author_cell_type=np.nan, braak_nft=str(m.braak_nft),
            path_braak_asyn_suppl=str(m.path_braak_asyn), cerad=str(m.cerad), gsm=r.gsm,
            ambient_removal="CellBender" if reg == "SN" else "DecontX")
        parts.append(a)
        print(r.file, a.shape, flush=True)
    a = ad.concat(parts, join="outer", fill_value=0)
    a.X = sp.csr_matrix(a.X)
    return a

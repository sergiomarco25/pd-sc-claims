"""Protein-level check of neuronal oxidative phosphorylation loss: imaging mass cytometry of nigral neurons.

Usage: python -I scripts/48_chen_imc_oxphos.py
Data (open, Figshare): Chen et al. 2021 (npj Parkinsons Dis; 10.25405/data.ncl.14079497), per-neuron log intensities of
OXPHOS complex I-V subunits (NDUFB8, NDUFA13, SDHA, UQCRC2, MTCO1, COXIV, ATP5B, OSCP) and mitochondrial mass markers
(TOM20, VDAC1) in substantia nigra neurons of PD and control cases; Chen et al. 2023 (10.25405/data.ncl.22005623),
per-neuron z-scores of mitochondrial proteins in TH-positive neurons (citrate retrieval). Per case, protein levels are
averaged over neurons (OXPHOS subunits also normalised to the mitochondrial mass marker VDAC1 / porin), and PD cases
are compared with controls (Mann-Whitney; Hedges' g on case means). Writes results/biology/chen_imc_oxphos.tsv.
"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

RAW = cfg.RAW
B = cfg.RESULTS / "biology"


def compare(cm, grp, case, ctrl, label):
    rows = []
    for c in cm.columns:
        a, r = cm.loc[grp == case, c].dropna(), cm.loc[grp == ctrl, c].dropna()
        sp = np.sqrt(((len(a) - 1) * a.var() + (len(r) - 1) * r.var()) / (len(a) + len(r) - 2))
        rows.append({"study": label, "protein": c, "n_pd": len(a), "n_ctrl": len(r),
                     "hedges_g": (a.mean() - r.mean()) / sp, "p": stats.mannwhitneyu(a, r).pvalue})
    return rows


def main():
    rows = []
    a = pd.read_csv(RAW / "chen2021" / "imc1_log.csv", encoding="utf-8-sig")
    a.columns = a.columns.str.strip()
    a["group"] = a.group.str.strip()
    prot = [c for c in a.columns if c.startswith("log_")]
    for c in prot:
        if c not in ("log_VDAC1", "log_TOM20"):
            a[c.replace("log_", "rel_")] = a[c] - a["log_VDAC1"]  # per VDAC1 (log scale)
    cm = a.groupby("caseno")[[c for c in a.columns if c.startswith(("log_", "rel_"))]].mean()
    grp = a.drop_duplicates("caseno").set_index("caseno").group.reindex(cm.index)
    rows += compare(cm, grp, "P", "Con", "Chen 2021, nigral neurons")
    b = pd.read_csv(RAW / "chen2023" / "citrate_z_TH_cat3.csv")
    b = b[b.grp.isin(["PD", "Control"])]
    z = [c for c in b.columns if c.startswith("z_") and not c.endswith(("DNA1", "Histone"))]
    cm = b.groupby("casecode")[z].mean()
    grp = b.drop_duplicates("casecode").set_index("casecode").grp.reindex(cm.index)
    rows += compare(cm, grp, "PD", "Control", "Chen 2023, TH+ neurons")
    t = pd.DataFrame(rows)
    t.to_csv(B / "chen_imc_oxphos.tsv", sep="\t", index=False)
    print(t.round(3).to_string(index=False))


if __name__ == "__main__":
    main()

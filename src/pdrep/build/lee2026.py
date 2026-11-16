"""lee2026: PsychAD (Lee et al. 2026 Nature), MSSM cohort CELLxGENE release, DLPFC only (4,140,453 nuclei, 1,042 donors).

Source: raw/psychad_mssm/0e853475-e298-4b09-881a-ed0b60d5a8c9.h5ad. X holds raw integer UMI counts (float32; no raw/
layers). Rows of kept donors are streamed in chunks.

Donor selection (lee2026 is used here for its Lewy body contrasts):
- keep donors whose CELLxGENE `disease` contains 'Parkinson disease' (Parkinson_disease flag = Yes) or
  'Lewy body dementia' (DLBD_status = Yes), whatever else they carry (AD, vascular, schizophrenia, PSP, MDD);
- keep controls = disease 'normal' (all have every *_status flag = No: neurotypical);
- drop everything else (AD-only, dementia NOS, schizophrenia, vascular, tauopathy, FTD, ...).
disease: 'Lewy body dementia' -> DLB (PsychAD DLBD); 'Parkinson disease' with 'dementia' -> PDD; 'Parkinson disease'
alone -> PD; 'normal' -> control. The two LB groups are disjoint in this release. `comorbidity` lists the other terms.
age: 'NN-year-old stage' -> NN; '80 year-old and over stage' -> 90 (CELLxGENE censoring; ages 81-89 appear explicitly,
so the censored donors are >= 90), flagged in `age_censored`. pmi / braak_lb: not in the release (NaN).
"""
import re
import numpy as np
import pandas as pd
import anndata as ad

from pdrep.harmonise import RAW
from pdrep.build._io import read_df, read_csr

PATH = RAW / "psychad_mssm" / "0e853475-e298-4b09-881a-ed0b60d5a8c9.h5ad"


def _cls(s):
    for pre, c in [("Oligo", "Oligodendrocytes"), ("OPC", "OPCs"), ("Astro", "Astrocytes"), ("Micro", "Microglia"),
                   ("EN_", "Excitatory neurons"), ("IN_", "Inhibitory neurons"), ("Endo", "Endothelial"),
                   ("PC_", "Pericytes"), ("SMC", "Pericytes"), ("VLMC", "Fibroblasts"), ("PVM", "Immune (other)"),
                   ("Adaptive_T", "T cells"), ("Adaptive", "Immune (other)")]:
        if s.startswith(pre):
            return c
    return "Unknown"


SUBTYPES = ["Oligo_OPALIN", "Oligo_RBFOX1", "Astro_WIF1", "Astro_PLSCR1", "Astro_ADAMTSL3", "Astro_GRIA1", "Micro",
            "OPC", "OPC_HAP1", "OPC_CHRM3", "OPC_GPR17", "EN_L2_3_IT_PDGFD", "EN_L2_3_IT_NTNG1", "EN_L3_5_IT_3",
            "EN_L3_5_IT_2_MET", "EN_L6_IT_1", "EN_L6_CT", "EN_L3_5_IT_1_PLSCR4", "EN_L3_5_IT_2_DACH1", "EN_L6B",
            "EN_L5_6_NP", "EN_L3_5_IT_1_CUX2", "EN_L6_IT_2", "EN_L6_CT_HSPA1A", "EN_L3_5_IT_3_HSPA1A",
            "EN_L6_IT_1_HSPA1A", "EN_L5_ET", "EN_L6_IT_2_HSPA1A", "EN_L3_5_IT_2_HSPA1A", "IN_VIP_TRPC6",
            "IN_PVALB_ANOS1", "IN_SST_EDNRA", "IN_VIP_BCL11B", "IN_LAMP5_RELN", "IN_PVALB_CHRM2", "IN_SST_EYA4",
            "IN_PVALB_CHC", "IN_LAMP5_LHX6", "IN_ADARB2_COL12A1", "IN_ADARB2_SV2C", "IN_ADARB2_SYT10", "IN_VIP_SCML4",
            "IN_SST_PRR16", "IN_ADARB2_RAB37", "IN_SST_MAML3", "IN_PVALB_HSPA1A", "IN_SST_NPY", "IN_LAMP5_RELN_HSPA1A",
            "IN_LAMP5_LHX6_HSPA1A", "IN_PVALB_STUM", "Endo_ITIH5", "Endo_IL1R1", "Endo_PCSK5", "PC_STAC", "PC_ADAMTS4",
            "SMC_MYOCD", "SMC_NRP1", "VLMC_ABCA6", "VLMC_DCDC2", "VLMC_SLC4A4", "PVM", "Adaptive_T", "Adaptive_NK",
            "Adaptive_B", "Adaptive_Plasma"]
AUTHOR_TO_CLASS = {s: _cls(s) for s in SUBTYPES}


def _disease(s):
    if s == "normal":
        return "control"
    if "Lewy body dementia" in s:
        return "DLB"
    if "Parkinson disease" in s:
        return "PDD" if "dementia" in s.split(" || ") else "PD"
    return None


def _age(s):
    m = re.match(r"^(\d+)-year-old", s)
    if m:
        return float(m.group(1))
    return 90.0 if s.startswith("80 year-old and over") else np.nan


def load():
    o = read_df(PATH, "obs")
    dis_raw = o["disease"].astype(str)
    dmap = {s: _disease(s) for s in dis_raw.unique()}
    disease = dis_raw.map(dmap)
    keep = disease.notna().values
    print(f"keeping {keep.sum()} / {len(keep)} cells, {o.donor_id[keep].nunique()} donors", flush=True)
    o, disease, dis_raw = o[keep], disease[keep], dis_raw[keep]
    var = read_df(PATH, "var")
    X = read_csr(PATH, "X", rows=keep, chunk=100_000)
    st = o["development_stage"].astype(str)
    lb_terms = {"dementia", "Parkinson disease", "Lewy body dementia"}
    donor = o["donor_id"].astype(str)
    obs = pd.DataFrame({
        "donor": donor.values, "sample": donor.values, "region": "DLPFC",
        "disease": disease.values, "disease_raw": dis_raw.values,
        "sex": o["sex"].astype(str).map({"male": "M", "female": "F"}).values,
        "age": st.map(_age).values, "age_censored": st.str.startswith("80 year-old and over").values,
        "pmi": np.nan, "braak_lb": np.nan,
        "author_cell_type": o["subtype"].astype(str).values, "subclass_author": o["subclass"].astype(str).values,
        "comorbidity": dis_raw.map(lambda s: " || ".join(t for t in s.split(" || ") if t not in lb_terms
                                                         and t != "normal") or "none").values,
        "ad_status": o["AD_status"].astype(str).values, "dlbd_status": o["DLBD_status"].astype(str).values,
        "parkinson_disease_flag": o["Parkinson_disease"].astype(str).values,
        "genetic_ancestry": o["genetic_ancestry"].astype(str).values, "brain_bank": "MSSM",
    }, index=o.index.astype(str))
    del o
    v = pd.DataFrame({"gene_ids": var.index.astype(str)}, index=pd.Index(var["feature_name"].astype(str).values))
    a = ad.AnnData(X, obs=obs, var=v)
    a.var_names_make_unique()
    return a

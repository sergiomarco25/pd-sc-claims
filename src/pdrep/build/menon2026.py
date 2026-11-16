"""menon2026: Menon et al. 2026 bioRxiv, Banner Sun Health multiome atlas; 101 donors (80 PD, 21 controls) x 5 regions
(SN, putamen, anterior cingulate, MTG, cerebellum), 4,502,568 nuclei (CELLxGENE collection 16876983-...).

Source: external_validation/menon2026/menon2026_complete_rna.h5ad. X is normalised; raw integer counts in raw/X
(same 25,183 genes), streamed in row chunks (logic copied from revision/scripts/33_menon_counts_cache.py).
Donor metadata (Supplementary Table 1): revision/results/30_menon2026/donor_metadata.csv (read-only).
- disease: Control=yes -> control; PD=yes and dementia_yn=1 -> PDD (spec: PDD only when the source says dementia);
  other PD -> PD. disease_raw = CELLxGENE disease label; `clinical_dx` = neurological_dx, `dementia` = dementia_yn.
- braak_lb: NaN. The cohort is staged with the Unified Staging System for LB disorders (not Braak), kept as
  `unified_lb_stage` (0, 2 = IIa, 3 = III, 4 = IV) and `unified_lb_stage_raw`; `braak_nft` from 'Braak score'.
"""
import numpy as np
import pandas as pd
import anndata as ad

from pdrep.harmonise import DATA
from pdrep.build._io import read_df, read_csr

PATH = DATA / "external_validation" / "menon2026" / "menon2026_complete_rna.h5ad"

from pdrep import config as cfg  # noqa: E402

META = str(cfg.LEGACY / "revision/results/30_menon2026/donor_metadata.csv")

AUTHOR_TO_CLASS = {
    "Oligo": "Oligodendrocytes", "OPC": "OPCs", "Astrocyte": "Astrocytes", "Bergmann-Glia": "Astrocytes",
    "Microglia": "Microglia", "PV-Macrophage": "Immune (other)", "T-cell": "T cells",
    "Excitatory": "Excitatory neurons", "RORB-Neuron": "Excitatory neurons",
    "Granule": "Excitatory neurons",  # cerebellar granule cells (glutamatergic)
    "DRD1-MSNs": "MSNs", "DRD2-MSNs": "MSNs", "DRD1-DRD2-MSNs": "MSNs",
    "Dopaminergic": "DA neurons",
    **{k: "Inhibitory neurons" for k in ["VIP", "PVALB", "SST", "LAMP5", "Chandelier", "GABA-Put", "GABA-SN",
                                           "MLI", "Golgi", "Purkinje"]},
    "Capillary": "Endothelial", "Pericyte": "Pericytes", "SMC": "Pericytes",
    "Fibroblast": "Fibroblasts", "PV-Fibroblast": "Fibroblasts", "Ependymal": "Ependymal",
}
REGION = {"SN": "SN", "PUT": "putamen", "CING": "cingulate", "MTG": "MTG", "CBL": "cerebellum"}
ROMAN = {"0": 0, "I": 1, "II": 2, "III": 3, "IV": 4, "V": 5, "VI": 6}


def _donor_meta():
    m = pd.read_csv(META).set_index("BannerID")
    ctrl, pd_, dem = m["Control"].eq("yes"), m["PD"].eq("yes"), m["dementia_yn"].eq(1)
    dis = pd.Series(np.select([ctrl, pd_ & dem, pd_], ["control", "PDD", "PD"], "other"), index=m.index)
    st = m["Unified LB Stage"].astype(str)
    uni = st.str.extract(r"^\s*(\w+)\.")[0].str.upper().map({"0": 0, "LLA": 2, "IIA": 2, "LLL": 3, "III": 3,
                                                              "LV": 4, "IV": 4})
    return pd.DataFrame({
        "disease": dis, "age": m["age_at_death"].astype(float), "pmi": m["PMI"].astype(float),
        "braak_lb": np.nan, "braak_nft": m["Braak score"].astype(str).str.strip().map(ROMAN),
        "unified_lb_stage": uni, "unified_lb_stage_raw": st, "clinical_dx": m["neurological_dx"].astype(str),
        "dementia": m["dementia_yn"].astype(float), "brain_bank": "Banner Sun Health",
    })


def load():
    meta = _donor_meta()
    assert meta.unified_lb_stage.notna().all(), meta.unified_lb_stage_raw.unique()
    o = read_df(PATH, "obs")
    var = read_df(PATH, "raw/var")
    X = read_csr(PATH, "raw/X", chunk=100_000)
    region = o["region"].astype(str).map(REGION)
    assert region.notna().all()
    donor = o["donor_id"].astype(str)
    obs = pd.DataFrame({
        "donor": donor.values, "sample": (donor + "_" + region).values, "region": region.values,
        "disease_raw": o["disease"].astype(str).values,
        "sex": o["sex"].astype(str).map({"male": "M", "female": "F"}).values,
        "author_cell_type": o["author_cell_type"].astype(str).values,
        "cell_type_cellxgene": o["cell_type"].astype(str).values,
        "batch": o.index.astype(str).str.split("#").str[0].values,
    }, index=o.index.astype(str))
    del o
    obs = obs.join(meta, on="donor")
    assert obs.disease.notna().all() and (obs.disease != "other").all()
    # consistency: CELLxGENE disease vs Suppl. Table 1
    assert ((obs.disease_raw == "normal") == (obs.disease == "control")).all()
    v = pd.DataFrame({"gene_ids": var.index.astype(str)}, index=pd.Index(var["feature_name"].astype(str).values))
    a = ad.AnnData(X, obs=obs, var=v)
    a.var_names_make_unique()
    return a

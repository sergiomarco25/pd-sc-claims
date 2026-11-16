"""prashant2024: AMP-PD multiregion snRNA-seq atlas (Sci Data 2024), 5 regions (DMNX, GPi, PMC, DLPFC, PVC).

Source: multiregion_sc/adata_multiregion_prashant.h5ad (CELLxGENE release, all 2,096,155 nuclei, 17,266 genes).
X is log-normalised; raw integer counts are in raw/X (same cells and genes; the *_2pct/_10pct[_raw] files are random
subsets of the same object and are not used). raw/X is streamed in row chunks to float32.
Author labels: obs.derived_class2 (12 classes incl. DMNX_Neu and GPI_Neu); CELLxGENE cell_type kept as extra column.
"""
import numpy as np
import pandas as pd
import anndata as ad

from pdrep.harmonise import DATA
from pdrep.build._io import read_df, read_csr

PATH = DATA / "multiregion_sc" / "adata_multiregion_prashant.h5ad"

AUTHOR_TO_CLASS = {
    "Oligo": "Oligodendrocytes", "EN": "Excitatory neurons", "Astro": "Astrocytes", "IN": "Inhibitory neurons",
    "Micro_PVM": "Microglia", "OPC": "OPCs", "Mural": "Pericytes", "Ependymal": "Ependymal",
    "Endo": "Endothelial", "Adaptive": "T cells",
    "DMNX_Neu": "Neurons (other)",  # DMNX (cholinergic/motor) neurons; CELLxGENE: 'central nervous system neuron'
    "GPI_Neu": "Neurons (other)",   # GPi projection neurons; CELLxGENE: 'central nervous system neuron'
}
TISSUE = {"prefrontal cortex": "DLPFC", "primary motor cortex": "motor", "primary visual cortex": "visual",
          "medial globus pallidus": "GPi", "dorsal motor nucleus of vagus nerve": "DMNX"}


def load():
    o = read_df(PATH, "obs")
    var = read_df(PATH, "raw/var")
    X = read_csr(PATH, "raw/X")
    region = o["tissue"].astype(str).map(TISSUE)
    assert region.notna().all()
    dis = o["disease"].astype(str)
    donor = o["donor_id"].astype(str)
    obs = pd.DataFrame({
        "donor": donor, "sample": donor + "_" + region, "region": region,
        "disease": dis.map({"Parkinson disease": "PD", "normal": "control"}), "disease_raw": dis,
        "sex": o["sex"].astype(str).map({"male": "M", "female": "F"}),
        "age": o["development_stage"].astype(str).str.extract(r"^(\d+)-year")[0].astype(float),
        "pmi": o["PMI"].astype(float), "braak_lb": o["path_braak_lb"].astype(float),
        "author_cell_type": o["derived_class2"].astype(str), "cell_type_cellxgene": o["cell_type"].astype(str),
        "brain_bank": o["Brain_bank"].astype(str), "rin": o["RIN"].astype(float),
        "ancestry": o["self_reported_ethnicity"].astype(str),
    })
    obs.index = o.index.astype(str)
    assert obs.disease.notna().all()
    v = pd.DataFrame({"gene_ids": var.index.astype(str)}, index=pd.Index(var["feature_name"].astype(str).values))
    a = ad.AnnData(X, obs=obs, var=v)
    a.var_names_make_unique()
    return a

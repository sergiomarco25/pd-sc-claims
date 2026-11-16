"""martirosyan2024: Martirosyan et al. 2024 Mol Neurodegener, SNpc snRNA-seq, 15 PD + 14 controls (Oregon Brain Bank).

Source: published_datasets/adata_martirosyan.h5ad (83,484 nuclei x 33,537 genes). X is a dense float32 matrix of raw
integer UMI counts (checked); read in row chunks and stored sparse. Author subpopulations in obs.SUBTYPES
(cell_type = broad class). One sample per donor (Sample ID s.0xxx, Brain Bank ID = donor).
Braak LB stage is not provided; Lewy body presence (midbrain / limbic / neocortical) and NFT Braak are kept.
"""
import h5py
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp

from pdrep.harmonise import DATA
from pdrep.build._io import read_df

PATH = DATA / "published_datasets" / "adata_martirosyan.h5ad"

AUTHOR_TO_CLASS = {
    **{f"oligos{i}": "Oligodendrocytes" for i in range(6)},
    **{f"astrocytes{i}": "Astrocytes" for i in range(6)},
    **{f"microglia{i}": "Microglia" for i in range(6)},
    "Astro": "Astrocytes", "Micro": "Microglia", "OPC": "OPCs",
    # Neurons0 = DA neurons (TH, SLC6A3, SLC18A2, ALDH1A1 markers, Suppl. Table neuron markers);
    # Neurons3 = GABAergic; others = non-DA neurons (see claims/signatures/martirosyan2024.tsv)
    "neurons00": "DA neurons", "neurons01": "Neurons (other)", "neurons02": "Neurons (other)",
    "neurons03": "Inhibitory neurons", "neurons04": "Neurons (other)", "neurons05": "Neurons (other)",
    "Neurons": "Neurons (other)",
    "VC": "Endothelial",  # vascular cells (endothelial + pericytes, not split by the authors)
    "T cells": "T cells",
}


def load():
    o = read_df(PATH, "obs")
    var = read_df(PATH, "var")
    with h5py.File(PATH, "r") as h:
        d = h["X"]
        blocks = []
        for s in range(0, d.shape[0], 5000):
            x = d[s:s + 5000]
            assert np.allclose(x[x > 0][:100000], np.round(x[x > 0][:100000]))
            blocks.append(sp.csr_matrix(x, dtype=np.float32))
    X = sp.vstack(blocks, format="csr")
    dx = o["Clinical diagnosis"].astype(str)
    obs = pd.DataFrame({
        "donor": o["Brain Bank ID"].astype(str), "sample": o["Sample ID"].astype(str), "region": "SN",
        "disease": dx.map({"Control": "control", "Parkinson's": "PD"}), "disease_raw": dx,
        "sex": o["Sex"].astype(str).map({"male": "M", "female": "F"}), "age": o["Age"].astype(float),
        "pmi": o["PMI hours"].astype(float), "braak_lb": np.nan,
        "author_cell_type": o["SUBTYPES"].astype(str), "author_broad": o["cell_type"].astype(str),
        "braak_nft": o["Braak stage for neurofibrillary tangles "].astype(float),
        "cerad": o["CERAD score for neuritic plaques"].astype(float),
        "lb_midbrain": o["Lewy bodies presence in midbrain"].astype(float),
        "lb_limbic": o["Lewy bodies presence in limbic regions (amygdala)"].astype(float),
        "lb_neocortex": o["Lewy bodies presence in neocortical regions (frontal cortex)"].astype(float),
        "rin": o["RIN measure"].astype(float), "brain_bank": "Oregon",
    })
    obs.index = o.index.astype(str)
    assert obs.disease.notna().all()
    return ad.AnnData(X, obs=obs, var=pd.DataFrame(index=var.index.astype(str)))

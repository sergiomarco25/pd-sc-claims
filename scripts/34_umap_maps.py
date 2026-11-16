"""Per-dataset UMAP maps with per-nucleus stress scores, to show the data behind the donor-level results.

Usage: python -I scripts/34_umap_maps.py <dataset>
Reads harmonised/<dataset>.h5ad (raw counts), takes up to 2,000 nuclei per donor x region and 40,000 in total
(seed 0; sorted nuclei excluded), normalises to 10,000 counts and log1p, selects 2,000 highly variable genes,
runs PCA (30 components), Harmony on donor, a 15-neighbour graph and UMAP. Per-nucleus scores (mean of the per-gene
z-scored log-normalised expression across the sampled nuclei): heat shock (HSPA1A, HSPA1B, HSPH1, DNAJB1, DNAJA1, HSP90AA1, BAG3, CHORDC1,
STIP1, AHSA1) and cold-inducible (CIRBP, RBM3); also log-normalised HSPA1A, CIRBP and SLC38A2.
Writes results/biology/maps/<dataset>.tsv.gz (coordinates, metadata, scores).
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import harmonypy

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

H = cfg.HARMONISED
OUT = cfg.RESULTS / "biology/maps"
HEAT = ["HSPA1A", "HSPA1B", "HSPH1", "DNAJB1", "DNAJA1", "HSP90AA1", "BAG3", "CHORDC1", "STIP1", "AHSA1"]
COLD = ["CIRBP", "RBM3"]
GENES = ["HSPA1A", "CIRBP", "SLC38A2"]
META = ["donor", "region", "disease", "sex", "age", "braak_lb", "cell_class"]


def main(ds):
    OUT.mkdir(parents=True, exist_ok=True)
    a = ad.read_h5ad(H / f"{ds}.h5ad", backed="r")
    obs = a.obs
    keep = np.ones(len(obs), bool)
    if "enrichment" in obs:
        keep &= ~obs.enrichment.astype(str).str.lower().eq("sorted").values
    rng = np.random.default_rng(0)
    pos = pd.Series(np.arange(len(obs)))[keep]
    grp = (obs.donor.astype(str) + "|" + obs.region.astype(str)).values[keep]
    take = pd.Series(pos.values, index=grp).groupby(level=0).apply(
        lambda x: x.sample(min(len(x), 2000), random_state=int(rng.integers(1e9))).values)
    idx = np.concatenate(take.values)
    if len(idx) > 40000:
        idx = rng.choice(idx, 40000, replace=False)
    idx = np.sort(idx)
    b = a[idx].to_memory()
    b.obs_names_make_unique()
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    try:
        sc.pp.highly_variable_genes(b, n_top_genes=2000, flavor="cell_ranger")
    except ValueError:  # cell_ranger binning fails on tied dispersion bins in some small datasets
        sc.pp.highly_variable_genes(b, n_top_genes=2000, flavor="seurat")
    sc.pp.pca(b, n_comps=30, mask_var="highly_variable")
    if b.obs.donor.nunique() > 1:
        ho = harmonypy.run_harmony(b.obsm["X_pca"], b.obs, "donor", max_iter_harmony=20, random_state=0)
        z = np.asarray(ho.Z_corr)
        b.obsm["X_pca_harmony"] = z if z.shape[0] == b.n_obs else z.T  # orientation differs between harmonypy versions
        rep = "X_pca_harmony"
    else:
        rep = "X_pca"
    sc.pp.neighbors(b, n_neighbors=15, use_rep=rep)
    sc.tl.umap(b, random_state=0)
    for name, gs in [("heat_score", HEAT), ("cold_score", COLD)]:  # mean of per-gene z-scored log expression
        g = [x for x in gs if x in b.var_names]
        v = b[:, g].X
        v = np.asarray(v.todense()) if hasattr(v, "todense") else np.asarray(v)
        b.obs[name] = ((v - v.mean(0)) / np.where(v.std(0) > 0, v.std(0), 1)).mean(1)
    out = b.obs[[m for m in META if m in b.obs]].copy()
    out["umap1"], out["umap2"] = b.obsm["X_umap"][:, 0], b.obsm["X_umap"][:, 1]
    out["heat_score"], out["cold_score"] = b.obs.heat_score.values, b.obs.cold_score.values
    for g in GENES:
        if g in b.var_names:
            v = b[:, g].X
            out[g] = np.asarray(v.todense()).ravel() if hasattr(v, "todense") else np.ravel(v)
    out["dataset"] = ds
    out.to_csv(OUT / f"{ds}.tsv.gz", sep="\t")
    print(ds, b.shape, out.cell_class.value_counts().to_dict())


if __name__ == "__main__":
    main(sys.argv[1])

"""Common schema for every harmonised dataset, and the steps all loaders share.

A loader (src/pdrep/build/<dataset>.py) defines `load() -> AnnData` with raw integer counts in X (cells x genes,
gene symbols as var_names) and the obs columns in OBS (values it cannot fill: NaN / "unknown"), plus a dict
`AUTHOR_TO_CLASS` mapping its `author_cell_type` values to CELL_CLASSES (empty when there are no author labels:
cells are then labelled by marker clustering here).

`finalise()` checks the schema, labels or checks cell classes against canonical markers, and writes
  harmonised/<dataset>.h5ad                 cells, raw counts, schema obs
  harmonised/<dataset>_donors.tsv           one row per donor x region (metadata + cells per class)
  harmonised/<dataset>_label_check.tsv      author label vs marker class agreement
  harmonised/pseudobulk/<dataset>.h5ad      summed counts per donor x region x cell_class (the unit for most tests)
"""
from pathlib import Path
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
import scanpy as sc

from pdrep import config as cfg  # noqa: E402

DATA = cfg.DATA_ROOT
RAW = DATA / "pd_sc_replication" / "raw"
OUT = DATA / "pd_sc_replication" / "harmonised"

OBS = ["dataset", "donor", "sample", "region", "disease", "disease_raw", "sex", "age", "pmi", "braak_lb",
       "author_cell_type", "cell_class"]
REGIONS = {"SN", "VTA", "midbrain", "putamen", "caudate", "striatum", "PFC", "DLPFC", "cingulate", "temporal", "MTG",
           "motor", "visual", "parietal", "amygdala", "hippocampus", "cerebellum", "DMNX", "GPi", "SVZ", "other"}
DISEASES = {"control", "PD", "PDD", "DLB", "LBD", "iLBD", "AD", "other"}
CELL_CLASSES = ["DA neurons", "MSNs", "Excitatory neurons", "Inhibitory neurons", "Neurons (other)", "Astrocytes",
                "Microglia", "Oligodendrocytes", "OPCs", "Endothelial", "Pericytes", "Fibroblasts", "T cells",
                "Immune (other)", "Ependymal", "Unknown"]
MARKERS = {
    "DA neurons": ["TH", "SLC6A3", "DDC", "SLC18A2", "KCNJ6"],
    "MSNs": ["PPP1R1B", "DRD1", "DRD2", "PENK", "TAC1"],
    "Excitatory neurons": ["SLC17A7", "SATB2", "SLC17A6", "CAMK2A", "NRGN"],
    "Inhibitory neurons": ["GAD1", "GAD2", "SLC32A1", "SST", "PVALB"],
    "Neurons (other)": ["RBFOX3", "SNAP25", "SYT1", "STMN2"],
    "Astrocytes": ["AQP4", "GJA1", "SLC1A2", "GFAP", "ALDH1L1", "SLC1A3"],
    "Microglia": ["P2RY12", "CX3CR1", "CSF1R", "DOCK8", "TMEM119"],
    "Oligodendrocytes": ["MBP", "PLP1", "MOBP", "ST18", "MOG"],
    "OPCs": ["PDGFRA", "VCAN", "PTPRZ1", "OLIG1"],
    "Endothelial": ["CLDN5", "FLT1", "VWF", "PECAM1"],
    "Pericytes": ["PDGFRB", "RGS5", "KCNJ8", "ACTA2"],
    "Fibroblasts": ["DCN", "COL1A2", "LUM", "COL1A1"],
    "T cells": ["CD3E", "CD3D", "CD2", "SKAP1", "CD247"],  # no THEMIS (L6 excitatory neurons) or IL7R
    "Immune (other)": ["MRC1", "F13A1", "CD163", "LYVE1", "MS4A1"],  # no PTPRC: microglia express it
    "Ependymal": ["FOXJ1", "RSPH1", "DYNLRB2", "CCDC153", "CFAP126"],  # not DNAH9/CFAP43: also in MSNs
}


def _norm_copy(a):
    b = a.copy()
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    return b


def marker_scores(a, max_cells=200_000, seed=0):
    """Per-cell score for each broad class (scanpy score_genes on log-normalised counts); subsampled if large."""
    if a.n_obs > max_cells:
        a = a[np.random.default_rng(seed).choice(a.n_obs, max_cells, replace=False)]
    b = _norm_copy(a)
    cols = {}
    for ct, genes in MARKERS.items():
        g = [x for x in genes if x in b.var_names]
        if len(g) >= 2:
            sc.tl.score_genes(b, g, score_name="s", random_state=seed, n_bins=25, ctrl_size=50)
            cols[ct] = b.obs["s"].values
    return pd.DataFrame(cols, index=b.obs_names)


NEURON_SUBTYPES = ["DA neurons", "MSNs", "Excitatory neurons", "Inhibitory neurons"]
LINEAGES = [c for c in MARKERS if c not in NEURON_SUBTYPES]  # "Neurons (other)" = pan-neuronal lineage score
# Classes only plausible in some regions (others may appear in any region).
REGION_ONLY = {"DA neurons": {"SN", "VTA", "midbrain"},
               "MSNs": {"putamen", "caudate", "striatum"},
               "Ependymal": {"SVZ", "caudate", "striatum"}}
SUBTYPE_MIN_Z = 0.5  # a neuron goes to a subtype only if that subtype's z-score clears this; else "Neurons (other)"
MIN_RAW = 0.1  # a cluster can only take a class whose raw mean marker score (above control genes) clears this
# Within the neuron step, these subtypes win over the listed competitor when they clear min_sub (MSNs are GABAergic).
PRECEDENCE = {"MSNs": "Inhibitory neurons"}
# A cluster can only take these classes if any listed gene is detected in at least this fraction of its nuclei
# (kamath2022 caudate MSNs scored "Ependymal" on DNAH9; wang2024 neuron/oligodendrocyte mixtures with TH in 12-15%
# of nuclei scored "DA neurons").
REQUIRE_DETECT = {"DA neurons": (["TH", "SLC6A3"], 0.30), "Ependymal": (["FOXJ1"], 0.20)}


def classify(z, region, min_sub=SUBTYPE_MIN_Z, raw=None):
    """Two-step class call from marker scores (rows = cells or clusters; z-scored for clusters, raw for cells).

    Step 1 picks the lineage (pan-neuronal vs each glial / vascular / immune class); step 2 picks the neuron subtype
    for neuron rows. Classes restricted by REGION_ONLY are excluded outside their regions. Z-scoring each class
    across cells first keeps one abundant class (e.g. ambient oligodendrocyte RNA in SN) from winning everywhere;
    `raw` (same shape, unscaled scores) then excludes classes whose markers are not actually expressed (raw < MIN_RAW),
    which stops near-absent classes (T cells, ependymal) from winning on noise.
    """
    if raw is not None:
        z = z.where(raw.reindex_like(z) >= MIN_RAW, -np.inf)
    region = pd.Series(np.asarray(region, dtype=object), index=z.index).astype(str)
    out = pd.Series("Unknown", index=z.index, dtype=object)
    for reg, idx in region.groupby(region).groups.items():
        zz = z.loc[idx]
        allowed = [c for c in zz.columns if c not in REGION_ONLY or reg in REGION_ONLY[c]]
        lin_scores = zz[[c for c in LINEAGES if c in allowed]]
        lin = lin_scores.idxmax(axis=1).where(np.isfinite(lin_scores.max(axis=1)), "Unknown")
        sub_cols = [c for c in NEURON_SUBTYPES if c in allowed]
        if sub_cols:
            best = zz[sub_cols].idxmax(axis=1)
            for win, lose in PRECEDENCE.items():
                if win in sub_cols:
                    best = best.where(~((best == lose) & (zz[win] >= min_sub)), win)
            ok = zz[sub_cols].max(axis=1) >= min_sub
            lin = lin.where(~((lin == "Neurons (other)") & ok), best)
        out.loc[idx] = lin.values
    return out


def zscore(s):
    return (s - s.mean()) / s.std(ddof=0).replace(0, 1)


def annotate_by_markers(a, resolution=1.0, seed=0):
    """For datasets without author labels: Leiden clusters labelled by `classify` on cluster-mean z-scored marker scores."""
    b = _norm_copy(a)
    sc.pp.highly_variable_genes(b, n_top_genes=3000)
    sc.pp.pca(b, n_comps=30, mask_var="highly_variable", random_state=seed)
    sc.pp.neighbors(b, n_neighbors=15, random_state=seed)
    sc.tl.leiden(b, resolution=resolution, key_added="leiden", flavor="igraph", n_iterations=2, random_state=seed)
    clusters = ("leiden_" + b.obs["leiden"].astype(str)).values
    return clusters, label_clusters(a, clusters)


def label_clusters(a, clusters):
    """Class per cluster from cluster-mean z-scored marker scores and the cluster's majority region."""
    s = marker_scores(a, max_cells=a.n_obs + 1)
    cl = pd.Series(np.asarray(clusters), index=a.obs_names)
    key = cl.loc[s.index].values
    zc, rc = zscore(s).groupby(key).mean(), s.groupby(key).mean()
    for cls, (genes, thr) in REQUIRE_DETECT.items():
        g = [x for x in genes if x in a.var_names]
        if cls in rc and g:
            det = pd.Series(np.asarray((a[s.index, g].X > 0).max(axis=1).todense()).ravel()
                            if hasattr(a.X, "todense") else (a[s.index, g].X > 0).max(axis=1), index=s.index)
            ok = det.groupby(key).mean() >= thr
            rc.loc[~ok.reindex(rc.index, fill_value=False), cls] = -np.inf
    reg = a.obs.region.astype(str).groupby(cl.values).agg(lambda r: r.value_counts().index[0])
    lab = classify(zc, reg.reindex(zc.index), raw=rc)
    return cl.map(lab).astype(str).values


def label_check(a):
    """For each author label: share of its cells whose marker class (`classify`, per cell) equals the mapped cell_class.

    `agree_lineage` counts any neuron class as agreeing with any other neuron class (neuron subtypes are hard to call
    per cell from five markers)."""
    s = marker_scores(a)  # raw per-cell scores: z-scoring single cells inflates noise in classes with no real cells
    obs = a.obs.loc[s.index]
    top = classify(s, obs.region, min_sub=0.0)
    neuron = set(NEURON_SUBTYPES) | {"Neurons (other)"}
    t = pd.DataFrame({"author_cell_type": obs.author_cell_type.astype(str).values,
                      "cell_class": obs.cell_class.astype(str).values, "marker_top": top.values})
    t["agree"] = t.marker_top == t.cell_class
    t["agree_lineage"] = t.agree | (t.marker_top.isin(neuron) & t.cell_class.isin(neuron))
    g = t.groupby(["author_cell_type", "cell_class"])
    return g.size().rename("n_cells_checked").to_frame().join(g[["agree", "agree_lineage"]].mean()).join(
        g.marker_top.agg(lambda x: x.value_counts().index[0]).rename("marker_majority")).reset_index()


def pseudobulk(a, keys=("donor", "region", "cell_class"), min_cells=10):
    """Sum raw counts per donor x region x cell_class (x enrichment, when sorted and unsorted nuclei are mixed)."""
    keys = list(keys) + (["enrichment"] if "enrichment" in a.obs and "enrichment" not in keys else [])
    grp = a.obs[list(keys)].astype(str).agg("|".join, axis=1)
    codes, uniq = pd.factorize(grp)
    ind = sp.csr_matrix((np.ones(a.n_obs), (codes, np.arange(a.n_obs))), shape=(len(uniq), a.n_obs))
    X = (ind @ sp.csr_matrix(a.X)).tocsr()
    meta_cols = [c for c in OBS if c not in ("author_cell_type", "sample")] + [k for k in keys if k not in OBS]
    first = a.obs.assign(_g=codes).groupby("_g")[meta_cols].first()
    obs = first.reindex(range(len(uniq)))
    obs["n_cells"] = np.bincount(codes)
    obs.index = uniq
    pb = ad.AnnData(X, obs=obs, var=a.var[[]].copy())
    return pb[pb.obs.n_cells >= min_cells].copy()


def finalise(a, dataset, author_to_class):
    """Validate, label, check and write a loader's output."""
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "pseudobulk").mkdir(exist_ok=True)
    a.obs["dataset"] = dataset
    for c in OBS:
        if c not in a.obs:
            a.obs[c] = np.nan
    X = a.X if sp.issparse(a.X) else sp.csr_matrix(a.X)
    probe = X[: min(5000, a.n_obs)].data
    assert probe.size and np.allclose(probe, np.round(probe)), f"{dataset}: X must hold raw integer counts"
    a.X = X.astype(np.float32, copy=False)  # no second copy of multi-GB matrices
    del X
    a.var_names_make_unique()
    bad = set(a.obs.region.dropna().astype(str)) - REGIONS
    assert not bad, f"{dataset}: regions not in vocabulary: {bad}"
    bad = set(a.obs.disease.dropna().astype(str)) - DISEASES
    assert not bad, f"{dataset}: diseases not in vocabulary: {bad}"

    if author_to_class:
        a.obs["cell_class"] = a.obs.author_cell_type.astype(str).map(author_to_class).fillna("Unknown")
        missing = set(a.obs.author_cell_type.astype(str)) - set(author_to_class)
        if missing:
            print(f"WARNING {dataset}: author labels without mapping -> Unknown: {sorted(missing)[:20]}")
    else:
        a.obs["author_cell_type"], a.obs["cell_class"] = annotate_by_markers(a)
    bad = set(a.obs.cell_class) - set(CELL_CLASSES)
    assert not bad, f"{dataset}: cell classes not in vocabulary: {bad}"
    sc.pp.calculate_qc_metrics(a, percent_top=None, log1p=False, inplace=True)

    for c in a.obs.columns:
        if a.obs[c].dtype == object:
            a.obs[c] = a.obs[c].astype(str).astype("category")
    a.write_h5ad(OUT / f"{dataset}.h5ad", compression="gzip")
    write_derived(a, dataset)
    return a


def write_derived(a, dataset):
    """Label check, donor table and pseudobulk for a finished cell-level object."""
    chk = label_check(a)
    chk.to_csv(OUT / f"{dataset}_label_check.tsv", sep="\t", index=False)
    unit = a.obs.groupby(["donor", "region"], observed=True)
    donors = unit[["disease", "disease_raw", "sex", "age", "pmi", "braak_lb"]].first().join(
        pd.crosstab([a.obs.donor, a.obs.region], a.obs.cell_class)).reset_index()
    donors.insert(0, "dataset", dataset)
    donors.to_csv(OUT / f"{dataset}_donors.tsv", sep="\t", index=False)
    pseudobulk(a).write_h5ad(OUT / "pseudobulk" / f"{dataset}.h5ad", compression="gzip")
    print(f"{dataset}: {a.n_obs} cells, {a.n_vars} genes, {a.obs.donor.nunique()} donors, "
          f"{len(donors)} donor x region units")
    print(donors.groupby(["region", "disease"]).donor.nunique().to_string())
    print(a.obs.cell_class.value_counts().to_string())
    w = chk.n_cells_checked / chk.n_cells_checked.sum()
    print("label agreement (cell-weighted): class", round((chk.agree * w).sum(), 3),
          "lineage", round((chk.agree_lineage * w).sum(), 3))


def relabel(dataset):
    """Re-run labelling checks on a built dataset after a change to the marker rules (no reloading or re-clustering).

    Datasets labelled by marker clustering (author_cell_type = leiden_*) get new cluster classes and a rewritten h5ad;
    author-labelled datasets keep their classes and only get a new label check, donor table and pseudobulk."""
    a = ad.read_h5ad(OUT / f"{dataset}.h5ad")
    if a.obs.author_cell_type.astype(str).str.startswith("leiden_").all():
        old = a.obs.cell_class.astype(str).copy()
        a.obs["cell_class"] = pd.Categorical(label_clusters(a, a.obs.author_cell_type.astype(str).values))
        print(f"{dataset}: relabelled clusters; old -> new class counts")
        print(pd.crosstab(old, a.obs.cell_class).to_string())
        a.write_h5ad(OUT / f"{dataset}.h5ad", compression="gzip")
    write_derived(a, dataset)

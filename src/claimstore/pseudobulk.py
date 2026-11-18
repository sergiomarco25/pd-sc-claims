"""Export one dataset's pseudobulk for edgeR: counts per region x cell class, and a unit table with one group column
per contrast (claimstore.contrasts, config/contrasts.toml).

Usage: python -I src/claimstore/cli.py pseudobulk DATASET      (analysis environment)
Writes results/de_cache/<dataset>/input/: <region>__<class>.csv.gz (genes x units, raw counts; class "ALL" = all
nuclei of the unit), units.tsv (unit metadata, donor-level extras of config/harmonise.toml, one column per contrast:
"case" / "ref" / empty), manifest.tsv. Datasets mixing sorted and unsorted nuclei (obs `enrichment`): the
`sorted_classes` of config/harmonise.toml use sorted units, every other class unsorted units. A dataset of sorted
nuclei only (e.g. a microglia release) has no "ALL": its nuclei are not the tissue.
"""
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

from claimstore import config as cfg
from claimstore.contrasts import CONTRASTS, unit_groups
from claimstore.harmonise import DONOR_EXTRAS, OUT as H, SORTED_CLASSES

RES = cfg.RESULTS / "de_cache"
MIN_CELLS = 10


def fname(cls):
    """File-name form of a cell class ("Neurons (other)" -> "Neurons_other")."""
    return cls.replace(" ", "_").replace("(", "").replace(")", "")


def donor_extras(ds):
    """Donor-level columns beyond the schema (config/harmonise.toml donor_extras) from the cell obs."""
    obs = ad.read_h5ad(H / f"{ds}.h5ad", backed="r").obs
    keep = [c for c in DONOR_EXTRAS if c in obs]
    return obs.groupby("donor", observed=True)[keep].first() if keep else pd.DataFrame(index=obs.donor.unique())


def export(ds):
    out = RES / ds / "input"
    out.mkdir(parents=True, exist_ok=True)
    pb = ad.read_h5ad(H / "pseudobulk" / f"{ds}.h5ad")
    pb.obs = pb.obs.join(donor_extras(ds), on="donor")
    tissue = True  # whether "ALL" stands for the tissue
    if "enrichment" in pb.obs:
        sorted_ = pb.obs.enrichment.astype(str).str.lower().eq("sorted")
        tissue = not sorted_.all()
        keep = np.where(pb.obs.cell_class.astype(str).isin(SORTED_CLASSES), sorted_, ~sorted_)
        pb = pb[keep].copy()
    pb.obs["unit"] = pb.obs.donor.astype(str) + "|" + pb.obs.region.astype(str)
    # ALL = sum over classes per unit (before the per-class n_cells filter)
    codes, units = pd.factorize(pb.obs.unit)
    ind = sp.csr_matrix((np.ones(pb.n_obs), (codes, np.arange(pb.n_obs))), shape=(len(units), pb.n_obs))
    allX = ind @ sp.csr_matrix(pb.X)
    meta = pb.obs.drop_duplicates("unit").set_index("unit").loc[units]
    meta["n_cells_ALL"] = np.bincount(codes, weights=pb.obs.n_cells.values)
    groups = unit_groups(meta)
    units_tab = meta.drop(columns=[c for c in ["cell_class", "n_cells"] if c in meta]).join(groups)
    units_tab.index.name = "unit"
    units_tab.to_csv(out / "units.tsv", sep="\t")

    rows = []
    genes = pb.var_names.values
    for reg in sorted(meta.region.astype(str).unique()):
        sel = meta.region.astype(str).values == reg
        blocks = [("ALL", allX[sel], units[sel], meta.n_cells_ALL.values[sel])] if tissue else []
        for cls in sorted(pb.obs.cell_class.astype(str).unique()):
            m = (pb.obs.region.astype(str).values == reg) & (pb.obs.cell_class.astype(str).values == cls) & \
                (pb.obs.n_cells.values >= MIN_CELLS)
            if m.sum() >= 4:
                blocks.append((cls, pb.X[np.where(m)[0]], pb.obs.unit.values[m], pb.obs.n_cells.values[m]))
        for cls, X, u, n in blocks:
            X = sp.csr_matrix(X)
            expressed = np.asarray((X > 0).sum(0)).ravel() >= 3
            df = pd.DataFrame(X[:, expressed].T.toarray().astype(np.int64), index=genes[expressed], columns=u)
            f = out / f"{reg}__{fname(cls)}.csv.gz"
            df.to_csv(f)
            ncell = pd.Series(n, index=u)
            for con in CONTRASTS:
                g = groups.loc[u, con]
                rows.append({"dataset": ds, "region": reg, "cell_class": cls, "file": f.name, "contrast": con,
                             "n_case": int((g == "case").sum()), "n_ref": int((g == "ref").sum()),
                             "median_cells": float(ncell.median())})
    man = pd.DataFrame(rows)
    man["run"] = (man.n_case >= 3) & (man.n_ref >= 3)
    man.to_csv(out / "manifest.tsv", sep="\t", index=False)
    print(ds, "fits to run:", int(man.run.sum()), "of", len(man))
    print(man[man.run].groupby(["contrast"]).size().to_string())


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) != 1 or argv[0].startswith("-"):
        print("usage: pseudobulk DATASET")
        return 2
    export(argv[0])
    return 0


if __name__ == "__main__":
    sys.exit(main())

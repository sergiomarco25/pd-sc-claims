"""Export one dataset's pseudobulk for edgeR: counts per region x cell class, and a unit table with one group column
per contrast (pdrep.contrasts.CONTRASTS).

Usage: python -I scripts/12_export_pseudobulk.py <dataset>
Writes results/de_cache/<dataset>/input/: <region>__<class>.csv.gz (genes x units, raw counts; class "ALL" = all
nuclei of the unit), units.tsv (unit metadata + one column per contrast: "case" / "ref" / empty), manifest.tsv.
Sorted-nuclei datasets (obs `enrichment`): DA neurons use sorted units, every other class unsorted units.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pdrep.harmonise import OUT as H  # noqa: E402
from pdrep.contrasts import CONTRASTS, unit_groups  # noqa: E402

from pdrep import config as cfg  # noqa: E402

RES = cfg.RESULTS / "de_cache"
MIN_CELLS = 10


def donor_extras(ds):
    """Donor-level columns beyond the schema (age_group, genotype, unified_lb_stage, ...) from the cell obs."""
    obs = ad.read_h5ad(H / f"{ds}.h5ad", backed="r").obs
    keep = [c for c in ["age_group", "genotype", "unified_lb_stage", "gliad_donor_group"] if c in obs]
    return obs.groupby("donor", observed=True)[keep].first() if keep else pd.DataFrame(index=obs.donor.unique())


def main(ds):
    out = RES / ds / "input"
    out.mkdir(parents=True, exist_ok=True)
    pb = ad.read_h5ad(H / "pseudobulk" / f"{ds}.h5ad")
    pb.obs = pb.obs.join(donor_extras(ds), on="donor")
    if "enrichment" in pb.obs:
        sorted_ = pb.obs.enrichment.astype(str).str.lower().eq("sorted")
        keep = np.where(pb.obs.cell_class.astype(str).eq("DA neurons"), sorted_, ~sorted_)
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
        blocks = [("ALL", allX[sel], units[sel], meta.n_cells_ALL.values[sel])]
        for cls in sorted(pb.obs.cell_class.astype(str).unique()):
            m = (pb.obs.region.astype(str).values == reg) & (pb.obs.cell_class.astype(str).values == cls) & \
                (pb.obs.n_cells.values >= MIN_CELLS)
            if m.sum() >= 4:
                blocks.append((cls, pb.X[np.where(m)[0]], pb.obs.unit.values[m], pb.obs.n_cells.values[m]))
        for cls, X, u, n in blocks:
            X = sp.csr_matrix(X)
            expressed = np.asarray((X > 0).sum(0)).ravel() >= 3
            df = pd.DataFrame(X[:, expressed].T.toarray().astype(np.int64), index=genes[expressed], columns=u)
            f = out / f"{reg}__{cls.replace(' ', '_').replace('(', '').replace(')', '')}.csv.gz"
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


if __name__ == "__main__":
    main(sys.argv[1])

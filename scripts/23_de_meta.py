"""Consensus gene expression changes in PD per cell type: random-effects meta-analysis of the pseudobulk edgeR fits.

Usage: python -I scripts/23_de_meta.py [contrast]   (default PDany_vs_control; also LBDany_vs_control, braak)
For each cell class, every dataset x region fit of the contrast (results/de_cache/*/de/) gives per gene a log2 fold
change and a standard error (|logFC| / |z|, z from the two-sided p). Regions of one dataset are combined by fixed
effect; datasets are pooled with DerSimonian-Laird (vectorised). Genes must be measured in at least half the datasets
of that class (minimum 3). Also per region group (midbrain / SN, striatum, cortex) where >= 3 datasets exist.
Writes results/biology/de_meta_<contrast>.tsv.gz (class, scope, gene, k, est, se, z, p, fdr, I2, n_pos, n_sig).
"""
from pathlib import Path
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.meta import bh  # noqa: E402
from pdrep.plot.style import REGION_GROUP  # noqa: E402

from pdrep import config as cfg  # noqa: E402

BASE = cfg.RESULTS
# --public: leave out the datasets of config [release] private_datasets and write to results/public/biology
PUBLIC = "--public" in sys.argv
DE, OUT = BASE / "de_cache", (BASE / "public" / "biology" if PUBLIC else BASE / "biology")


def fname_to_class(stem):
    return stem.replace("_", " ").replace("Neurons other", "Neurons (other)").replace("Immune other", "Immune (other)")


def load(contrast):
    rows = []
    for ds_dir in sorted(DE.iterdir()):
        log = ds_dir / "de_log.tsv"
        if not log.exists() or os.path.getsize(log) < 60 or (PUBLIC and not cfg.is_public(ds_dir.name)):
            continue
        lg = pd.read_csv(log, sep="\t")
        for r in lg[lg.contrast == contrast].itertuples():
            f = ds_dir / "de" / f"{r.file[:-len('.csv.gz')]}__{contrast}.tsv.gz"
            if not f.exists():
                continue
            t = pd.read_csv(f, sep="\t", usecols=["gene", "logFC", "p"])
            z = stats.norm.isf(np.clip(t.p.values, 1e-300, 1) / 2)
            se = np.abs(t.logFC.values) / np.maximum(z, 1e-3)
            rows.append(pd.DataFrame({"dataset": ds_dir.name, "region": r.region, "cell_class": r.cell_class,
                                      "gene": t.gene.str.upper().values, "est": t.logFC.values, "se": se,
                                      "sig": t.p.values < 0.05, "n": r.n_case + (r.n_ref if pd.notna(r.n_ref) else 0)}))
    return pd.concat(rows, ignore_index=True)


def dl_matrix(E, S):
    """Vectorised DerSimonian-Laird over rows (genes) of est/se matrices with NaNs."""
    W = 1 / S ** 2
    ok = np.isfinite(E) & np.isfinite(W)
    E0, W0 = np.where(ok, E, 0), np.where(ok, W, 0)
    k = ok.sum(1)
    sw = W0.sum(1)
    mu_f = (W0 * E0).sum(1) / np.where(sw > 0, sw, np.nan)
    q = (W0 * (E0 - mu_f[:, None]) ** 2 * ok).sum(1)
    c = sw - (W0 ** 2).sum(1) / np.where(sw > 0, sw, np.nan)
    tau2 = np.where((k > 1) & (c > 0), np.maximum(0, (q - (k - 1)) / np.where(c > 0, c, 1)), 0)
    WR = np.where(ok, 1 / (S ** 2 + tau2[:, None]), 0)
    mu = (WR * E0).sum(1) / np.where(WR.sum(1) > 0, WR.sum(1), np.nan)
    se = np.sqrt(1 / np.where(WR.sum(1) > 0, WR.sum(1), np.nan))
    z = mu / se
    i2 = np.where((k > 1) & (q > 0), np.maximum(0, (q - (k - 1)) / np.where(q > 0, q, 1)), 0)
    return mu, se, z, 2 * stats.norm.sf(np.abs(z)), i2, k


def meta_block(d, min_k):
    # combine regions within dataset (fixed effect), then pool datasets
    d = d.assign(w=1 / d.se ** 2)
    g = d.groupby(["gene", "dataset"])
    per = pd.DataFrame({"est": g.apply(lambda x: (x.w * x.est).sum() / x.w.sum()) if False else
                        (d.w * d.est).groupby([d.gene, d.dataset]).sum() / g.w.sum(),
                        "se": np.sqrt(1 / g.w.sum()), "sig": g.sig.any()})
    E = per.est.unstack()
    S = per.se.unstack().reindex_like(E)
    keep = E.notna().sum(1) >= min_k
    E, S = E[keep], S[keep]
    mu, se, z, p, i2, k = dl_matrix(E.values, S.values)
    sig = per.sig.unstack().reindex_like(E).fillna(False)
    return pd.DataFrame({"gene": E.index, "k": k, "est": mu, "se": se, "z": z, "p": p, "I2": i2,
                         "n_pos": (E.values > 0).sum(1), "n_sig_same_dir": ((sig.values) & (np.sign(E.values) == np.sign(mu)[:, None])).sum(1)})


def main(contrast):
    OUT.mkdir(parents=True, exist_ok=True)
    d = load(contrast)
    d = d[np.isfinite(d.est) & np.isfinite(d.se) & (d.se > 0)]
    d["region_group"] = d.region.map(REGION_GROUP).fillna("other")
    res = []
    for cls, dc in d.groupby("cell_class"):
        n_ds = dc.dataset.nunique()
        if n_ds >= 3:
            r = meta_block(dc, max(3, int(np.ceil(n_ds / 2))))
            res.append(r.assign(cell_class=cls, scope="all regions", n_datasets=n_ds))
        for rg, dr in dc.groupby("region_group"):
            n = dr.dataset.nunique()
            if n >= 3 and rg != "other":
                r = meta_block(dr, max(3, int(np.ceil(n / 2))))
                res.append(r.assign(cell_class=cls, scope=rg, n_datasets=n))
        print(cls, n_ds, "datasets", flush=True)
    out = pd.concat(res, ignore_index=True)
    out["fdr"] = out.groupby(["cell_class", "scope"]).p.transform(lambda p: pd.Series(bh(p.values), index=p.index))
    out.to_csv(OUT / f"de_meta_{contrast}.tsv.gz", sep="\t", index=False)
    s = out[out.scope == "all regions"].groupby("cell_class").apply(
        lambda x: pd.Series({"datasets": x.n_datasets.iloc[0], "genes": len(x), "fdr05_up": int(((x.fdr < 0.05) & (x.est > 0)).sum()),
                             "fdr05_down": int(((x.fdr < 0.05) & (x.est < 0)).sum())}))
    print(s.to_string())


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--public"]
    main(args[0] if args else "PDany_vs_control")

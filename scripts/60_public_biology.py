"""Public biology tables for the public atlas (claims/COMMUNITY_PLAN.md step 5.1).

Usage: python -I scripts/60_public_biology.py [--check]
Rebuilds every cross-dataset table the atlas shows without the datasets of config [release] private_datasets and
without the private papers' gene modules ([release] private_modules), into results/public/biology and
results/public/psychad. Each table is re-pooled from the per-dataset results with the pooling code of the script that
made it, so it is what that script gives on the public datasets alone:
  composition_meta.tsv        22 pool() on composition_estimates.tsv
  programme_meta.tsv          27 pool() on programme_effects.tsv (FDR over the public programmes)
  genetics_meta.tsv           25 pool() on genetics_per_dataset.tsv
  donor_heterogeneity.tsv     30 heterogeneity() on donor_scores.tsv.gz
  donor_regions_effects.tsv   30 regions() on the same; region_map.tsv by 39 pool_sites()
  psychad/specificity.tsv     41b rows without the private modules, FDR recomputed (PsychAD itself is public)
  consensus_genes_PDany.tsv   28 consensus() on the public gene meta-analysis, once 23_de_meta.py --public has run
The gene meta-analysis and the pathways run separately: 23_de_meta.py --public, 24_pathways.py --public.
--check first re-pools the full per-dataset tables and requires the stored full tables back, to show that the
pooling functions are the ones that made them.
"""
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep import config as cfg  # noqa: E402

B, P = cfg.RESULTS / "biology", cfg.RESULTS / "psychad"
OUT, OUTP = cfg.RESULTS / "public" / "biology", cfg.RESULTS / "public" / "psychad"


def script(name):
    spec = importlib.util.spec_from_file_location("s" + name[:3].rstrip("_"), ROOT / "scripts" / name)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def read(path, **kw):
    return pd.read_csv(path, sep="\t", low_memory=False, **kw)


def public(t):
    keep = t.dataset.map(cfg.is_public)
    if "programme" in t:
        keep &= ~t.programme.isin(cfg.PRIVATE_MODULES)
    return t[keep]


def same(name, got, want, keys):
    """Exact agreement (to 1e-9) of a recomputed table with the stored one, matched on keys."""
    got = got.sort_values(keys).reset_index(drop=True)
    want = want.sort_values(keys).reset_index(drop=True)
    assert len(got) == len(want), f"{name}: {len(got)} rows, stored {len(want)}"
    assert list(got[keys].astype(str).values.ravel()) == list(want[keys].astype(str).values.ravel()), f"{name}: keys differ"
    for c in want.columns:
        if c not in got:
            raise AssertionError(f"{name}: column {c} missing")
        a, b = got[c], want[c]
        if pd.api.types.is_numeric_dtype(b) and pd.api.types.is_numeric_dtype(a):
            ok = np.isclose(a.astype(float), b.astype(float), rtol=1e-9, atol=1e-12, equal_nan=True)
        else:
            ok = a.fillna("").astype(str).values == b.fillna("").astype(str).values
        assert ok.all(), f"{name}: column {c} differs in {int((~ok).sum())} rows"
    print(f"check {name}: {len(want)} rows identical")


def tables(s22, s25, s27, s28, s30, s39, flt):
    """All public (flt=public) or full (flt=identity) tables from the per-dataset inputs."""
    est = flt(read(B / "composition_estimates.tsv"))
    eff = flt(read(B / "programme_effects.tsv"))
    per = flt(read(B / "genetics_per_dataset.tsv"))
    sc = flt(read(B / "donor_scores.tsv.gz"))
    sc = sc[sc.programme.isin(s30.CORE)]
    het_ds, het = s30.heterogeneity(sc)
    _, _, reg_eff = s30.regions(sc)
    out = {"composition_estimates.tsv": est, "composition_meta.tsv": s22.pool(est),
           "programme_effects.tsv": eff, "programme_meta.tsv": s27.pool(eff),
           "genetics_per_dataset.tsv": per, "genetics_meta.tsv": s25.pool(per),
           "donor_heterogeneity.tsv": het, "donor_heterogeneity_per_dataset.tsv": het_ds,
           "donor_regions_effects.tsv": reg_eff, "region_map.tsv": s39.pool_sites(reg_eff)}
    spec = read(P / "specificity.tsv")
    if flt is public:
        spec = spec[~spec.programme.isin(cfg.PRIVATE_MODULES)].copy()
    spec["fdr"] = stats.false_discovery_control(spec.p)
    out["specificity.tsv"] = spec
    return out


KEYS = {"composition_meta.tsv": ["contrast", "cell_class", "denominator", "scope"],
        "programme_meta.tsv": ["programme", "cell_class", "test"], "genetics_meta.tsv": ["set_type", "cell_class"],
        "donor_heterogeneity.tsv": ["programme", "cell_class"], "region_map.tsv": ["programme", "cell_class", "site"],
        "donor_regions_effects.tsv": ["dataset", "programme", "cell_class", "region"],
        "specificity.tsv": ["programme", "cell_class", "group"]}


def main():
    s22, s25, s27 = script("22_composition_meta.py"), script("25_genetics_meta.py"), script("27_programmes.py")
    s28, s30, s39 = script("28_biology_figures.py"), script("30_donor_analyses.py"), script("39_region_map.py")
    if "--check" in sys.argv:
        full = tables(s22, s25, s27, s28, s30, s39, lambda t: t)
        for name, keys in KEYS.items():
            same(name, full[name], read((P if name == "specificity.tsv" else B) / name), keys)
        same("consensus_genes_PDany.tsv", s28.consensus(read(B / "de_meta_PDany_vs_control.tsv.gz")),
             read(B / "consensus_genes_PDany.tsv"), ["gene", "cell_class"])
    OUT.mkdir(parents=True, exist_ok=True)
    OUTP.mkdir(parents=True, exist_ok=True)
    pub = tables(s22, s25, s27, s28, s30, s39, public)
    for name, t in pub.items():
        assert not t.get("dataset", pd.Series(dtype=str)).isin(cfg.PRIVATE_DATASETS).any(), name
        assert not t.get("programme", pd.Series(dtype=str)).isin(cfg.PRIVATE_MODULES).any(), name
        t.to_csv((OUTP if name == "specificity.tsv" else OUT) / name, sep="\t", index=False)
        print(f"public {name}: {len(t)} rows")
    dm = OUT / "de_meta_PDany_vs_control.tsv.gz"
    if dm.exists():
        cons = s28.consensus(read(dm))
        cons.to_csv(OUT / "consensus_genes_PDany.tsv", sep="\t", index=False)
        print(f"public consensus_genes_PDany.tsv: {len(cons)} rows")
    else:
        print("consensus genes wait for: 23_de_meta.py PDany_vs_control --public")


if __name__ == "__main__":
    main()

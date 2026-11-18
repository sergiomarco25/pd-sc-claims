"""Cross-dataset meta-analysis of every claim: combine its tests in all other datasets (claimstore engine).

Usage: python -I src/claimstore/cli.py meta
Reads results/tests/cross/*.tsv (and internal/*.tsv for the source result); writes results/meta/claim_meta.tsv
(one row per claim x variant x scope) and results/meta/claim_dataset_matrix.tsv (claim x dataset outcomes).

Orientation: gene-level and state tests already report z with positive = as claimed; composition and concordance
tests ([tests] concordance in config/project.toml) report the raw sign, so they are multiplied by the claim's expected
sign. Combined z: weighted Stouffer
(weights sqrt(units)); heterogeneity: Q = sum (z_i - mean z)^2 with k-1 df, I2 = max(0, (Q - (k-1)) / Q).
Scopes: all other datasets; same region only; other regions only; datasets published after the claim's paper.
Robustness class (scope = all):
  robust          k >= 3, combined p < 0.05 as claimed, and >= 50% of datasets reproduce it
  supported       combined p < 0.05 as claimed (fewer than half reproduce it individually, or k = 2)
  contradicted    combined p < 0.05 in the opposite direction
  inconsistent    combined p >= 0.05 and I2 >= 50% (datasets disagree)
  not_replicated  combined p >= 0.05 and I2 < 50%
  insufficient    tested in fewer than 2 other datasets
"""
import glob
import sys

import numpy as np
import pandas as pd
import yaml
from scipy import stats

from claimstore import config as cfg

ROOT = cfg.REPO
RES = cfg.RESULTS
UP, DOWN = {"up", "enriched", "associated_pos", "exists"}, {"down", "depleted", "associated_neg"}
RAW_SIGN_TEMPLATES = {"composition_clr"} | {v["template"] for v in cfg.CONFIG.get("tests", {}).get("concordance", {}).values()}


def load(mode):
    """Tests of `mode` plus the separate GWAS pass (<mode>_gwas), which replaces any GWAS rows of the main pass."""
    parts = []
    for m in (mode, mode + "_gwas"):
        fs = [f for f in glob.glob(str(RES / "tests" / m / "*.tsv")) if open(f).read().count("\n") > 1]
        if fs:
            parts.append(pd.concat([pd.read_csv(f, sep="\t") for f in fs], ignore_index=True).assign(_part=m))
    if not parts:
        return pd.DataFrame(columns=["claim_id", "dataset", "outcome", "variant"])
    d = pd.concat(parts, ignore_index=True)
    d["variant"] = d["variant"].fillna("both") if "variant" in d else "both"
    gw = d[d._part.str.endswith("_gwas")]
    d = d[~(d._part.eq(mode) & d.set_index(["claim_id", "dataset"]).index.isin(gw.set_index(["claim_id", "dataset"]).index))]
    return d.drop(columns="_part")


def combine(g):
    t = g[g.outcome != "not_testable"].dropna(subset=["z_oriented"])
    k = len(t)
    out = {"k": k, **{f"n_{o}": int((t.outcome == o).sum()) for o in
                      ["reproduced", "direction_only", "not_reproduced", "contradicted"]}}
    if k == 0:
        return out | {"z": np.nan, "p": np.nan, "I2": np.nan, "frac_reproduced": np.nan, "frac_direction": np.nan}
    w = np.sqrt(t.n_units.clip(lower=1).values)
    z = t.z_oriented.values
    zc = (w * z).sum() / np.sqrt((w ** 2).sum())
    q = ((z - z.mean()) ** 2).sum()
    i2 = max(0.0, (q - (k - 1)) / q) if k > 1 and q > 0 else 0.0
    return out | {"z": zc, "p": 2 * stats.norm.sf(abs(zc)), "I2": i2,
                  "frac_reproduced": out["n_reproduced"] / k,
                  "frac_direction": (out["n_reproduced"] + out["n_direction_only"]) / k,
                  "datasets": ";".join(sorted(t.dataset))}


def robustness(r, null_claim=False):
    if r["k"] < 2:
        return "insufficient"
    if null_claim:  # "no difference" claims: equivalence-style; per-dataset outcomes already mapped by the engine
        if r["p"] < 0.05:
            return "contradicted"
        return "robust" if r["k"] >= 3 and r["frac_reproduced"] >= 0.5 else "supported" if r["frac_reproduced"] >= 0.5 \
            else "not_replicated"
    if r["p"] < 0.05 and r["z"] > 0:
        return "robust" if r["k"] >= 3 and r["frac_reproduced"] >= 0.5 else "supported"
    if r["p"] < 0.05 and r["z"] < 0:
        return "contradicted"
    return "inconsistent" if r["I2"] >= 0.5 else "not_replicated"


def main():
    claims = pd.DataFrame(yaml.safe_load((ROOT / "claims" / "claims.yaml").read_text())).set_index("id")
    papers = pd.read_csv(ROOT / "inventory" / "papers.tsv", sep="\t", dtype=str).set_index("paper_id")
    ctx = pd.read_csv(ROOT / "claims" / "contexts.tsv", sep="\t").set_index("claim_id")
    nulls = ctx["null_claim"].astype(str).str.lower().eq("true").to_dict() if "null_claim" in ctx else {}
    x = load("cross")
    if x.empty:
        print("no cross-dataset tests yet (results/tests/cross/)")
        return 0
    for c in ("z", "p", "effect", "n_case", "n_ref"):  # absent when no claim was testable
        if c not in x:
            x[c] = np.nan
    if "region_match" not in x:
        x["region_match"] = ""
    x["direction"] = x.claim_id.map(claims.direction)
    sign = x.direction.map(lambda v: 1 if v in UP else -1 if v in DOWN else 1)
    x["z_oriented"] = np.where(x.template_used.isin(RAW_SIGN_TEMPLATES), x.z * sign, x.z)
    x["n_units"] = x.n_case.fillna(0) + x.n_ref.fillna(0)
    x["claim_year"] = x.claim_id.map(ctx.paper).map(papers.year).astype(float)
    x["dataset_year"] = x.dataset.map(lambda d: papers.year.get(d)).astype(float)  # a dataset id is its paper's id

    rows = []
    for (cid, var), g in x.groupby(["claim_id", "variant"]):
        for scope, sub in [("all", g), ("same_region", g[g.region_match == "same"]),
                           ("other_region", g[g.region_match == "other"]),
                           ("later_datasets", g[g.dataset_year > g.claim_year])]:
            r = combine(sub)
            r.update({"claim_id": cid, "variant": var, "scope": scope})
            r["robustness"] = robustness(r, bool(nulls.get(cid, False))) if scope == "all" else ""
            rows.append(r)
    meta = pd.DataFrame(rows)
    out = RES / "meta"
    out.mkdir(exist_ok=True)
    meta.to_csv(out / "claim_meta.tsv", sep="\t", index=False)
    mat = x.pivot_table(index=["claim_id", "variant"], columns="dataset", values="outcome", aggfunc="first")
    mat.to_csv(out / "claim_dataset_matrix.tsv", sep="\t")
    allv = meta[(meta.scope == "all") & (meta.variant != "revised")]
    print(allv.robustness.value_counts().to_string())
    print(meta.groupby("scope").k.describe()[["count", "mean", "50%", "max"]].to_string())


if __name__ == "__main__":
    sys.exit(main())

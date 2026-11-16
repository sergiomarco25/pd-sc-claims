"""Build reports/R4_cross.html: replication of every claim across the other datasets.

Usage: python -I scripts/17_build_r4_report.py   (after scripts/15_meta_analysis.py)
Reads results/meta/claim_meta.tsv, results/meta/claim_dataset_matrix.tsv, internal results (via 15_meta_analysis.load),
claims/claims.yaml, inventory/papers.tsv. Headline view: GLIAD preprint variant ('revised' shown separately).
"""
from pathlib import Path
import importlib.util
import json
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

RES = cfg.RESULTS
CLASSES = ["robust", "supported", "inconsistent", "not_replicated", "contradicted", "insufficient"]
spec = importlib.util.spec_from_file_location("meta", ROOT / "scripts" / "15_meta_analysis.py")
meta_mod = importlib.util.module_from_spec(spec)
spec.loader.exec_module(meta_mod)


def unit_class(u):
    u = str(u).lower()
    return "cells" if u.startswith("cells") or u == "rois" else "donors / samples" if ("donor" in u or "sample" in u) else "other"


def pick_published(df):
    """One row per claim: GLIAD preprint variant first, then 'both'."""
    df = df[df.variant != "revised"]
    return df.sort_values("variant", key=lambda s: s.map({"preprint": 0, "both": 1})).drop_duplicates("claim_id")


def main():
    claims = pd.DataFrame(yaml.safe_load((ROOT / "claims" / "claims.yaml").read_text())).set_index("id")
    papers = pd.read_csv(ROOT / "inventory" / "papers.tsv", sep="\t", dtype=str).set_index("paper_id")
    ctx = pd.read_csv(ROOT / "claims" / "contexts.tsv", sep="\t").set_index("claim_id")
    m = pd.read_csv(RES / "meta" / "claim_meta.tsv", sep="\t")
    allm = pick_published(m[m.scope == "all"]).set_index("claim_id")
    scope = {s: pick_published(m[m.scope == s]).set_index("claim_id") for s in ["same_region", "other_region", "later_datasets"]}
    internal = pick_published(meta_mod.load("internal")).set_index("claim_id").outcome

    t = allm.copy()
    t["paper"] = t.index.map(ctx.paper)
    t["year"] = t.paper.map(papers.year)
    for col in ["statement", "centrality", "claim_type", "direction", "unit_of_analysis", "method", "quote", "location"]:
        t[col] = t.index.map(claims[col])
    t["unit_class"] = t.unit_of_analysis.map(unit_class)
    t["template"] = t.index.map(ctx.test_template)
    t["comparison"] = t.index.map(ctx.comparison_text)
    t["internal"] = t.index.map(internal).fillna("not_testable")
    for s, d in scope.items():
        t[f"{s}_k"] = t.index.map(d.k).fillna(0).astype(int)
        t[f"{s}_z"] = t.index.map(d.z)
        t[f"{s}_p"] = t.index.map(d.p)
        t[f"{s}_rep"] = t.index.map(d.frac_reproduced)
    tested = t[t.robustness != "insufficient"]

    def by(key, frame=tested):
        x = pd.crosstab(frame[key], frame.robustness).reindex(columns=CLASSES[:5], fill_value=0)
        x["n"] = x.sum(1)
        return [{"key": str(k), **{c: int(v) for c, v in r.items()}} for k, r in x.sort_values("n", ascending=False).iterrows()]

    rate = lambda col: [{"scope": s, "claims": int((t[f"{s}_k"] > 0).sum()),
                         "mean_frac_reproduced": float(np.nanmean(t.loc[t[f"{s}_k"] > 0, f"{s}_rep"])) if (t[f"{s}_k"] > 0).any() else None,
                         "combined_sig": int(((t[f"{s}_p"] < 0.05) & (t[f"{s}_z"] > 0)).sum())}
                        for s in ["same_region", "other_region", "later_datasets"]]
    cross_int = pd.crosstab(tested.internal, tested.robustness).reindex(columns=CLASSES[:5], fill_value=0)

    mat = pd.read_csv(RES / "meta" / "claim_dataset_matrix.tsv", sep="\t")
    mat = mat[mat.variant != "revised"].sort_values("variant", key=lambda s: s.map({"preprint": 0, "both": 1})).drop_duplicates("claim_id")
    mat = mat.set_index("claim_id").drop(columns="variant")
    mat = mat.loc[[i for i in tested.sort_values(["robustness", "paper"]).index if i in mat.index]]
    gl = m[(m.scope == "all") & m.claim_id.str.startswith("barbareyes2026")].pivot_table(
        index="claim_id", columns="variant", values="robustness", aggfunc="first").reindex(columns=["preprint", "revised", "both"]).fillna("")
    gl["statement"] = gl.index.map(claims.statement)
    gw = t[t.template.eq("genetic_enrichment")]

    cols = ["paper", "year", "statement", "centrality", "claim_type", "direction", "comparison", "template", "internal",
            "robustness", "k", "n_reproduced", "n_direction_only", "n_not_reproduced", "n_contradicted", "z", "p", "I2",
            "frac_reproduced", "datasets", "same_region_k", "same_region_rep", "other_region_k", "other_region_rep",
            "later_datasets_k", "later_datasets_rep", "unit_of_analysis", "method", "quote", "location", "unit_class"]
    rows = t.reindex(columns=cols).reset_index().rename(columns={"index": "claim_id"})
    data = {
        "claims": json.loads(rows.to_json(orient="records")),
        "counts": {c: int((t.robustness == c).sum()) for c in CLASSES} | {"n": len(t)},
        "by_paper": by("paper"), "by_template": by("template"), "by_unit": by("unit_class"),
        "by_centrality": by("centrality"), "by_type": by("claim_type"),
        "scopes": rate(None),
        "cross_internal": [{"internal": i, **{c: int(v) for c, v in r.items()}} for i, r in cross_int.iterrows()],
        "matrix": {"datasets": list(mat.columns), "rows": [{"claim_id": i, "robustness": t.loc[i, "robustness"],
                                                             "cells": [None if pd.isna(v) else v for v in r.values]}
                                                            for i, r in mat.iterrows()]},
        "gliad": json.loads(gl.reset_index().to_json(orient="records")),
        "gwas": json.loads(gw.reset_index().rename(columns={"index": "claim_id"})[
            ["claim_id", "statement", "robustness", "k", "n_reproduced", "internal", "datasets"]].to_json(orient="records")),
    }
    tpl = (ROOT / "scripts" / "r4_template.html").read_text()
    out = ROOT / "reports" / "R4_cross.html"
    out.write_text(tpl.replace("/*DATA*/null", json.dumps(data, ensure_ascii=False)))
    rows.to_csv(ROOT / "reports" / "R4_cross_claims.tsv", sep="\t", index=False)
    print("wrote", out, data["counts"])


if __name__ == "__main__":
    main()

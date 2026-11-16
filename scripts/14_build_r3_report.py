"""Build reports/R3_internal.html: internal replication of every claim in its own dataset.

Usage: python -I scripts/14_build_r3_report.py
Reads results/tests/internal/*.tsv, claims/claims.yaml, claims/contexts.tsv, claims/source_discrepancies.tsv.
Headline numbers use the as-published view (GLIAD preprint variant; 'revised' rows shown separately).
"""
from pathlib import Path
import glob
import json
import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

RES = cfg.RESULTS / "tests/internal"
OUTCOMES = ["reproduced", "direction_only", "not_reproduced", "contradicted", "not_testable"]


def unit_class(u):
    u = str(u).lower()
    return "cells" if u.startswith("cells") or u == "rois" else "donors / samples" if ("donor" in u or "sample" in u) else "other"


def why_untestable(note):
    n = str(note)
    if "later pass" in n:
        return "claim type not yet covered (communication, GWAS, trajectory, chromatin, spatial)"
    if "no signature genes" in n:
        return "no gene list available (supplement missing or text-only)"
    if "enough units" in n or "too few nuclei" in n:
        return "too few donors, nuclei or genes in the source data"
    if n in ("not_testable",):
        return "not testable by design (e.g. morphology, IHC, cross-modality)"
    if n == "not_available":
        return "needs data no dataset has (spatial tier, white matter, haplotype, IHC score)"
    return "other"


def main():
    import importlib.util
    spec = importlib.util.spec_from_file_location("meta", ROOT / "scripts" / "15_meta_analysis.py")
    meta = importlib.util.module_from_spec(spec); spec.loader.exec_module(meta)
    d = meta.load("internal")
    claims = pd.DataFrame(yaml.safe_load((ROOT / "claims" / "claims.yaml").read_text())).set_index("id")
    ctx = pd.read_csv(ROOT / "claims" / "contexts.tsv", sep="\t").set_index("claim_id")
    disc = pd.read_csv(ROOT / "claims" / "source_discrepancies.tsv", sep="\t")
    papers = pd.read_csv(ROOT / "inventory" / "papers.tsv", sep="\t", dtype=str).set_index("paper_id")

    pub = d[d.variant != "revised"].copy()
    # one row per claim for the as-published view: GLIAD preprint row wins over 'both'
    pub = pub.sort_values("variant", key=lambda s: s.map({"preprint": 0, "both": 1})).drop_duplicates("claim_id")
    pub["paper"] = pub.claim_id.str.rsplit("-C", n=1).str[0]
    for col, src in [("statement", "statement"), ("centrality", "centrality"), ("claim_type", "claim_type"),
                     ("quote", "quote"), ("location", "location"), ("authors_unit", "unit_of_analysis"),
                     ("authors_method", "method"), ("reported_stat", "reported_stat"), ("direction", "direction")]:
        pub[col] = pub.claim_id.map(claims[src])
    pub["unit_class"] = pub.authors_unit.map(unit_class)
    pub["comparison"] = pub.claim_id.map(ctx.comparison_text)
    pub["why_untestable"] = np.where(pub.outcome == "not_testable", pub.note.map(why_untestable), "")
    pub["year"] = pub.paper.map(papers.year)
    tested = pub[pub.outcome != "not_testable"]

    def table(key):
        x = pd.crosstab(tested[key], tested.outcome).reindex(columns=OUTCOMES[:4], fill_value=0)
        x["n"] = x.sum(1)
        return [{"key": k, **{o: int(v) for o, v in r.items()}} for k, r in x.iterrows()]

    by_paper = pd.crosstab(pub.paper, pub.outcome).reindex(columns=OUTCOMES, fill_value=0)
    gl = d[d.claim_id.str.startswith("barbareyes2026")].pivot_table(index="claim_id", columns="variant", values="outcome",
                                                                   aggfunc="first")
    gl = gl.reindex(columns=["preprint", "revised", "both"])
    gl = gl[gl[["preprint", "revised"]].notna().any(axis=1)].fillna("")
    gl["statement"] = gl.index.map(claims.statement)
    disc["statement"] = disc.claim_id.map(claims.statement).fillna("(applies to the whole paper)")

    cols = ["claim_id", "paper", "year", "statement", "centrality", "claim_type", "direction", "comparison", "outcome",
            "template_used", "z", "p", "effect", "n_case", "n_ref", "n_genes", "concordance", "frac_sig_concordant",
            "donor_spread", "regions_tested", "classes_tested", "authors_unit", "authors_method", "reported_stat",
            "quote", "location", "note", "why_untestable", "unit_class"]
    rows = pub.reindex(columns=cols)
    for c in ["z", "p", "effect", "concordance", "frac_sig_concordant", "donor_spread"]:
        rows[c] = pd.to_numeric(rows[c], errors="coerce").round(4)
    data = {
        "claims": json.loads(rows.to_json(orient="records")),
        "counts": {o: int((pub.outcome == o).sum()) for o in OUTCOMES} | {"n": len(pub), "papers": int(pub.paper.nunique())},
        "by_paper": [{"paper": p, **{o: int(v) for o, v in r.items()}} for p, r in by_paper.iterrows()],
        "by_template": table("template_used"), "by_unit": table("unit_class"), "by_centrality": table("centrality"),
        "untestable": pub[pub.outcome == "not_testable"].why_untestable.value_counts().to_dict(),
        "gliad": json.loads(gl.reset_index().to_json(orient="records")),
        "discrepancies": json.loads(disc.to_json(orient="records")),
    }
    tpl = (ROOT / "scripts" / "r3_template.html").read_text()
    out = ROOT / "reports" / "R3_internal.html"
    out.write_text(tpl.replace("/*DATA*/null", json.dumps(data, ensure_ascii=False)))
    print("wrote", out, data["counts"])
    pub.to_csv(ROOT / "reports" / "R3_internal_claims.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()

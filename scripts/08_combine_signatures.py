"""Combine per-paper signatures into claims/signatures_all.tsv.gz and per-claim coverage into claims/signature_coverage.tsv.

Usage: python -I scripts/08_combine_signatures.py

Gene symbols stay as the authors reported them in `gene_reported`; `gene` holds a cleaned symbol:
- Excel date corruption reversed: "2022-03-01 ..." / Excel serial 44256 / "1-Mar" → MARCH1; month 9 → SEPT<d>, month 12 → DEC<d>
  (legacy symbols; mapping to current HGNC names, e.g. MARCHF1 / SEPTIN7, happens per dataset at test time);
- non-breaking spaces and trailing ".00" removed ("BOP\xa01.00" → BOP1).
Coverage: every claim needing a signature (state_signature, geneset_shift, de_pseudobulk, or composition with a cell
state) is classed full (in a signature, no unresolved row), partial (both), text_only (only named_genes or
text-derived signatures), or none (unresolved only).
"""
from pathlib import Path
import datetime as dt
import re
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]
SIG = ROOT / "claims" / "signatures"
MONTH = {3: "MARCH", 9: "SEPT", 12: "DEC"}


def clean_gene(g):
    s = str(g).replace("\xa0", "").strip()
    d = None
    m = re.fullmatch(r"(\d{4})-(\d{2})-(\d{2})(?: 00:00:00)?", s)
    if m:
        d = dt.date(int(m[1]), int(m[2]), int(m[3]))
    elif re.fullmatch(r"4\d{4}", s):  # Excel serial day number (1900 system)
        d = dt.date(1899, 12, 30) + dt.timedelta(days=int(s))
    m2 = re.fullmatch(r"(\d{1,2})-(Mar|Sep|Dec)", s, re.I)  # Excel display form, e.g. "1-Mar"
    if m2:
        return f"{MONTH[{'mar': 3, 'sep': 9, 'dec': 12}[m2[2].lower()]]}{int(m2[1])}"
    if d is not None and d.month in MONTH:
        return f"{MONTH[d.month]}{d.day}"
    return re.sub(r"\.00$", "", s)


def main():
    frames, unres = [], []
    for f in sorted(SIG.glob("*.tsv")):
        if f.stem.endswith(("_unresolved", "_regions")):
            continue
        df = pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False)
        if "variant" not in df:
            df["variant"] = "both"
        frames.append(df)
        u = SIG / f"{f.stem}_unresolved.tsv"
        if u.exists():
            unres.append(pd.read_csv(u, sep="\t", dtype=str, keep_default_na=False).assign(paper=f.stem))
    sig = pd.concat(frames, ignore_index=True)
    sig["gene_reported"] = sig.gene
    sig["gene"] = sig.gene_reported.map(clean_gene)
    fixed = sig[sig.gene != sig.gene_reported]
    sig.to_csv(ROOT / "claims" / "signatures_all.tsv.gz", sep="\t", index=False)

    claims = pd.DataFrame(yaml.safe_load((ROOT / "claims" / "claims.yaml").read_text())).fillna("")
    need = claims[claims.test_template.isin(["state_signature", "geneset_shift", "de_pseudobulk"]) |
                  ((claims.claim_type == "composition") & (claims.cell_state != ""))]
    per = sig.assign(claim_id=sig.claim_ids.str.split(";")).explode("claim_id")
    per["claim_id"] = per.claim_id.str.strip()
    text_based = per.kind.eq("named_genes") | per.source.str.lower().eq("text") | per.selection.str.contains("text", case=False)
    agg = per.groupby("claim_id").agg(n_signatures=("signature_id", "nunique"),
                                      n_genes=("gene", "nunique"),
                                      table_based=("signature_id", lambda s: (~text_based[s.index]).any()))
    u = pd.concat(unres, ignore_index=True)
    reasons = u.groupby("claim_id").reason.agg(lambda r: " | ".join(sorted(set(r))))
    cov = need[["id", "paper", "claim_type", "test_template", "centrality"]].rename(columns={"id": "claim_id"})
    cov = cov.merge(agg, left_on="claim_id", right_index=True, how="left").fillna({"n_signatures": 0, "n_genes": 0})
    cov["unresolved_reason"] = cov.claim_id.map(reasons).fillna("")
    cov["coverage"] = [("none" if n == 0 else "text_only" if not tb else "partial" if r else "full")
                       for n, tb, r in zip(cov.n_signatures, cov.table_based.astype(object).where(cov.table_based.notna(), False), cov.unresolved_reason)]
    cov.to_csv(ROOT / "claims" / "signature_coverage.tsv", sep="\t", index=False)

    print(f"{sig.signature_id.nunique()} signatures, {len(sig)} gene rows, {sig.paper.nunique()} papers")
    print(f"gene symbols repaired: {fixed.gene_reported.nunique()} distinct, e.g.",
          dict(list(zip(fixed.gene_reported, fixed.gene))[:6]))
    print("\nclaims needing a signature:", len(cov))
    print(cov.coverage.value_counts().to_string())
    print("\n" + pd.crosstab(cov.paper, cov.coverage).to_string())


if __name__ == "__main__":
    main()

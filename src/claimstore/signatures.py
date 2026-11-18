"""Claim gene signatures for the re-tests (claimstore.retest): claims/signatures_all.tsv.gz, one row per claim x gene.

Two sources, combined:
  curated   claims/signatures/<paper>.tsv files a project writes from the papers' supplementary tables (columns at
            least signature_id, paper, claim_ids, kind, cell_type, gene; kind de_up / de_down / state_markers /
            named_genes / pathway_genes ...; claim_ids ';'-separated). Files ending in _unresolved.tsv are skipped.
  named     the genes each claim names in its `genes` field (claims/claims.yaml), kind named_genes, for every claim
            without a curated row. A named gene's expected direction is the claim's own.
A curated signature replaces the named genes of the claims it covers. Projects with their own combining step (the
reference project's scripts/08_combine_signatures.py) do not need this command.

Usage: python3 -I src/claimstore/cli.py signatures
"""
import sys

import pandas as pd
import yaml

from claimstore import config as cfg

ROOT = cfg.REPO
OUT = ROOT / "claims" / "signatures_all.tsv.gz"
# the columns the store reads (claimstore.build_store), in the reference project's order
COLUMNS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison", "gene", "rank", "stat",
           "padj", "source", "selection", "variant", "gene_reported"]


def as_list(v):
    if v is None or (isinstance(v, float) and v != v):
        return []
    if isinstance(v, (list, tuple)):
        return [str(x).strip() for x in v if str(x).strip()]
    return [x.strip() for x in str(v).replace(";", ",").split(",") if x.strip()]


def named(claims, covered=()):
    """named_genes rows from the claims' `genes` fields, for claims not in `covered`."""
    rows = []
    for c in claims:
        if c["id"] in covered:
            continue
        for g in dict.fromkeys(as_list(c.get("genes"))):  # unique, in order
            rows.append({"signature_id": f"{c['id']}_named", "paper": c.get("paper", ""), "claim_ids": c["id"],
                         "kind": "named_genes", "cell_type": ";".join(as_list(c.get("cell_type"))),
                         "cell_state": str(c.get("cell_state") or ""), "comparison": str(c.get("comparison") or ""),
                         "gene": g, "rank": "", "stat": "", "padj": "", "source": "claim genes field",
                         "selection": "named in the claim text", "variant": "", "gene_reported": g})
    return pd.DataFrame(rows, columns=COLUMNS)


def curated(folder):
    """Every claims/signatures/<paper>.tsv, with the columns of COLUMNS (missing ones empty)."""
    fs = sorted(f for f in folder.glob("*.tsv") if not f.name.endswith("_unresolved.tsv")) if folder.is_dir() else []
    if not fs:
        return pd.DataFrame(columns=COLUMNS)
    d = pd.concat([pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False) for f in fs], ignore_index=True)
    for c in COLUMNS:
        if c not in d:
            d[c] = d["gene"] if c == "gene_reported" else ""
    return d[COLUMNS]


def build(claims, cur):
    covered = {i.strip() for ids in cur.claim_ids for i in str(ids).split(";") if i.strip()}
    return pd.concat([cur, named(claims, covered)], ignore_index=True)


def main(argv=None):
    claims = yaml.safe_load((ROOT / "claims" / "claims.yaml").read_text()) or []
    cur = curated(ROOT / "claims" / "signatures")
    sig = build(claims, cur)
    sig.to_csv(OUT, sep="\t", index=False, compression="gzip")
    n_claims = sig.claim_ids.str.split(";").explode().str.strip().nunique() if len(sig) else 0
    print(f"{OUT.relative_to(ROOT)}: {len(sig)} rows, {sig.signature_id.nunique()} signatures, {n_claims} claims "
          f"({len(cur)} curated rows, {int(sig.kind.eq('named_genes').sum() - cur.kind.eq('named_genes').sum())} "
          "named genes from the claims)")
    return 0


if __name__ == "__main__":
    sys.exit(main())

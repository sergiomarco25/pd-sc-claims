"""Run claim tests in one dataset.

Usage: python -I scripts/13_run_tests.py internal <dataset>        claims whose source dataset is <dataset>
       python -I scripts/13_run_tests.py cross <dataset>           all other papers' claims, in <dataset>
       add "gwas" as a third argument to run only GWAS-enrichment claims (written to <mode>_gwas/), so they can
       run once the dataset's MAGMA gene sets exist without holding up everything else
Writes results/tests/<mode>/<dataset>.tsv (without GWAS claims when run separately; meta-analysis merges both).
"""
from pathlib import Path
import sys
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.tests.engine import test_claim, RES, H  # noqa: E402
from pdrep.regions import target_regions  # noqa: E402
import yaml  # noqa: E402

if __name__ == "__main__":
    mode, ds = sys.argv[1], sys.argv[2]
    part = sys.argv[3] if len(sys.argv) > 3 else "all"
    ctx = pd.read_csv(ROOT / "claims" / "contexts.tsv", sep="\t").fillna("")
    sig = pd.read_csv(ROOT / "claims" / "signatures_all.tsv.gz", sep="\t", dtype=str, keep_default_na=False)
    sig = sig.assign(claim_id=sig.claim_ids.str.split(";")).explode("claim_id")
    sig["claim_id"] = sig.claim_id.str.strip()
    sig["gene"] = sig.gene.where(sig.gene != "", sig.gene_reported)
    ctx = ctx[ctx.test_template.eq("genetic_enrichment")] if part == "gwas" else \
        ctx[~ctx.test_template.eq("genetic_enrichment")] if part == "nogwas" else ctx
    claims_all = pd.DataFrame(yaml.safe_load((ROOT / "claims" / "claims.yaml").read_text())).set_index("id")
    todo = ctx[ctx.dataset.eq(ds)] if mode == "internal" else ctx[~ctx.dataset.eq(ds)]
    if mode == "cross":  # retarget regions to this dataset; skip claims that are not testable anywhere
        claims = pd.DataFrame(yaml.safe_load((ROOT / "claims" / "claims.yaml").read_text())).set_index("id")
        have = set(pd.read_csv(H / f"{ds}_donors.tsv", sep="\t").region.astype(str))
        todo = todo[todo.status.eq("ok")].copy()
        reg = todo.claim_id.map(lambda i: target_regions(claims.region.get(i) if isinstance(claims.region.get(i), list)
                                                       else [claims.region.get(i)], have))
        todo["regions"], todo["region_match"] = reg.str[0].str.join(";"), reg.str[1]
    rows = []
    for c in todo.itertuples():
        if c.status != "ok":
            rows.append({"claim_id": c.claim_id, "dataset": ds, "outcome": "not_testable", "note": c.status})
            continue
        extra = {"region_match": c.region_match} if mode == "cross" else {}
        for r in test_claim(c, ds, sig[sig.claim_id.eq(c.claim_id)], claims_all.loc[c.claim_id]):
            rows.append({"claim_id": c.claim_id, "dataset": ds, "contrast": c.contrast, **extra, **r})
    out = RES / "tests" / (mode + ("_gwas" if part == "gwas" else ""))
    out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows, columns=None if rows else ["claim_id", "dataset", "outcome"])
    df.to_csv(out / f"{ds}.tsv", sep="\t", index=False)
    print(ds, mode, len(df), "rows")
    if len(df):
        print(df.outcome.value_counts().to_string())

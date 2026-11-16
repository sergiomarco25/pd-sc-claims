"""Gene signatures for prashant2024 (AMP-PD multiregion snRNA-seq + WGS data descriptor, Sci Data 2024).

Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/prashant2024.py
The paper has one claim (C01: fewer nuclei per sample in early-affected regions), a sample-yield claim with no genes,
cell state or pathway. No signature is needed: the script writes an empty claims/signatures/prashant2024.tsv (header
only) and records the claim in claims/signatures/prashant2024_unresolved.tsv.
"""
from pathlib import Path
import pandas as pd
import yaml

PAPER = "prashant2024"
ROOT = Path(__file__).resolve().parents[2]
OUTDIR = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison", "gene", "rank", "stat",
        "padj", "source", "selection"]

sig = pd.DataFrame(columns=COLS)
un = pd.DataFrame([dict(claim_id=f"{PAPER}-C01", term="nuclei yield per sample by region",
                        reason="no gene signature needed: claim is about nuclei counts per sample (metadata only); "
                               "not_testable template, no cell state or genes")])
OUTDIR.mkdir(parents=True, exist_ok=True)
sig.to_csv(OUTDIR / f"{PAPER}.tsv", sep="\t", index=False)
un.to_csv(OUTDIR / f"{PAPER}_unresolved.tsv", sep="\t", index=False)

claims = yaml.safe_load(open(ROOT / "claims" / f"{PAPER}.yaml"))
need = [c["id"] for c in claims if c["test_template"] in ("state_signature", "geneset_shift", "de_pseudobulk")
        or ((c["test_template"] == "composition_donor" or c["claim_type"] == "composition") and (c.get("cell_state") or "").strip())]
missing = [c for c in need if c not in set(un.claim_id)]
print(f"{len(claims)} claims; 0 signatures; {len(need)} claims need a signature; recorded in unresolved: "
      f"{list(un.claim_id)}; missing: {missing or 'none'}")

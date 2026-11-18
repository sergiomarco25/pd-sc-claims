"""Write the claim schemas (COMMUNITY_PLAN steps 2.3 and 2.6) from config/schema.toml with the engine's builder.

Outputs JSON Schema (draft-07) files next to this script:
  claim.schema.json     one published claim in claims/<paper_id>.yaml (claims/SCHEMA.md describes each field)
  finding.schema.json   one finding of this study in claims/findings.yaml (same fields plus origin and evidence)

Validated in CI by tests/test_claim_schema.py. Change config/schema.toml, run this script, commit both files.
Usage: python3 -I claims/schema/build_schemas.py      (standard library only)
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "src"))
from claimstore.schema import *  # noqa: E402,F401,F403  (SCHEMA_VERSION, VOCAB, REQUIRED, TEXT, claim_schema, ...)
from claimstore.schema import write  # noqa: E402

if __name__ == "__main__":
    for p in write(HERE):
        print(p)

"""Command line for a claim-store project: one entry point for the engine's steps.

Usage: python3 -I src/claimstore/cli.py COMMAND [ARGS...]
  schema      write claims/schema/*.json from config/schema.toml
  merge       validate the per-paper claim files and merge them (claims/claims.yaml, claims.tsv, validation.tsv)
  contexts    map every claim to its dataset, regions, cell classes and contrast (claims/contexts.tsv)
  supplements [--list] [PAPER...]: list (--list) or download the papers' supplementary files (raw/supplements/)
  signatures  claim gene signatures for the re-tests: curated claims/signatures/*.tsv plus the genes each claim names
  harmonise   DATASET: run the dataset's loader, check the schema, label cell classes, write harmonised files
  relabel     DATASET...: re-run labelling checks, donor table and pseudobulk after a change to the rules
  pseudobulk  DATASET: export counts per region x cell class and the unit table for edgeR
  de-cache    DATASET...: edgeR fits per region x class x contrast and for the stage (R: [envs] r_env)
  de-log      DATASET...: rebuild de_log.tsv from written fit files
  tests       internal|cross DATASET [gwas|nogwas]: re-test claims in one dataset (results/tests/)
  meta        combine every claim's tests across datasets (results/meta/claim_meta.tsv)
  store       build the claim store (results/store/claim_store.sqlite, claims/store/)
  rating      rate every claim of the store (claims/store/claim_rating.tsv); --db / --out
  atlas       the browsable atlas page of the store (reports/atlas.html); --public --db STORE for a release
  release     build a public release (VERSION ...)
  serve       MCP server over stdio;  http: MCP + REST server over HTTP
  config      print paths (--show, --env, NAME)
The steps that need pandas (contexts onwards) run in the project's analysis environment; harmonise, relabel and
pseudobulk also need anndata and scanpy, de-cache needs R with edgeR and data.table.
"""
import importlib
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

COMMANDS = {"schema": "schema", "merge": "claims", "contexts": "contexts", "signatures": "signatures", "supplements": "supplements", "harmonise": "harmonise",
            "relabel": "harmonise:relabel_main", "pseudobulk": "pseudobulk", "de-cache": "de_cache",
            "de-log": "de_cache:rebuild_main", "tests": "retest", "meta": "claim_meta", "store": "build_store",
            "rating": "rating", "atlas": "atlas", "release": "release", "serve": "mcp_server", "http": "mcp_http", "config": "config",
            "diff": "release_diff", "linked-data": "linked_data", "audit": "audit"}


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] in ("-h", "--help") or argv[0] not in COMMANDS:
        print(__doc__)
        return 0 if argv and argv[0] in ("-h", "--help") else 2
    name, _, fn = COMMANDS[argv[0]].partition(":")
    mod = importlib.import_module(f"claimstore.{name}")
    sys.argv = [f"claimstore {argv[0]}", *argv[1:]]
    out = getattr(mod, fn or "main")()
    return out if isinstance(out, int) else 0


if __name__ == "__main__":
    sys.exit(main())

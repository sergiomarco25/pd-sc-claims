# PD single-cell claims: a re-tested evidence base

Published conclusions of human post-mortem brain Parkinson's disease single-cell and single-nucleus studies, extracted
as structured claims, each re-tested donor-level in its own data and in every other harmonised dataset, and rated for
certainty (verdict and GRADE-style certainty; rules in `claims/RATING.md`). The same store holds this study's own
findings, rated the same way, consensus statements and reference knowledge, and is served as an atlas, a knowledge
graph and a server for AI tools (Model Context Protocol).

**Status:** release 0.9.1-preview (428 claims, 18 datasets, built 2026-11-16).

## What is here
| Path | What |
|---|---|
| `release/` | the current public release: `claim_store.sqlite`, every table as TSV, the knowledge graph (GraphML, JSON), claims as schema.org JSON-LD, an RO-Crate description, `release.json` and `SHA256SUMS` |
| `docs/index.html` | the atlas, a single static page built from the release, served at https://sergiomarco25.github.io/pd-sc-claims/ |
| `claims/` | the claims of each published paper (`<paper>.yaml`, schema in `claims/schema/`), signatures, contexts, consensus statements, reference knowledge and the rating rules |
| `src/pdrep/`, `scripts/`, `workflow/` | the pipeline: harmonising datasets, donor-level tests, meta-analysis, rating, the store, the atlas and the release (Snakemake; `config/project.toml` holds every path) |
| `deploy/mcp/` | the MCP server: run it locally against `release/claim_store.sqlite` |
| `tests/` | tests on a small public test store (run by CI) |
| `template/`, `scripts/63_new_project.py`, `.claude/skills/claim-replication/` | a claim store for another field: the project template, the command that makes a project from it, and the method step by step as a Claude Code skill |

One harmonised cohort whose own paper is not yet published is not part of this release: its data, its results and its
paper's claims are left out, and every rating is computed without it. It will be added after publication.

## Use
- Browse: the atlas at https://sergiomarco25.github.io/pd-sc-claims/, or open `docs/index.html`.
- Query: `sqlite3 release/claim_store.sqlite`, or the TSV files in `release/tables/`.
- AI tools: see `deploy/mcp/` to serve the store over MCP.
- Rebuild: the pipeline needs the raw datasets (accessions in `inventory/`); `snakemake -n` shows the plan.
- Ask a question of the evidence base: open an issue with the "Ask the evidence base" form.
- Another field: `python3 -I scripts/63_new_project.py DIR --name NAME --title TITLE` makes a new project; in Claude
  Code, the `claim-replication` skill (`.claude/skills/`) walks through the steps from papers to a public release.

## Licence and citation
Code: Apache-2.0 (`LICENSE`). Claims, store, tables and other data and text created here: CC BY 4.0 (`LICENSE-DATA.md`).
Harmonised data keep their sources' licences. Cite with `CITATION.cff`.

# Changelog

## 2026-11-18 — The engine, the template for another field and the claim-replication skill
- `src/claimstore/`: the field-independent engine (configuration, contrasts, harmonisation, pseudobulk, edgeR DE cache,
  claim re-tests, meta-analysis, rating, store, knowledge graph, atlas, release, linked data, MCP servers,
  supplementary files and curated signatures); `src/pdrep/` keeps the Parkinson's instance. Facts about the field are
  configuration: `config/schema.toml`, `contrasts.toml`, `contexts.toml`, `harmonise.toml` and `project.toml`.
- `template/` and `scripts/63_new_project.py`: a claim store for another field in one command, with its own checks
  and a copy of the engine.
- `.claude/skills/claim-replication/`: the method step by step as a Claude Code skill, from papers to a public release,
  with the claim extraction brief; every new project gets a copy.
- The release (0.9.1-preview) and the atlas are unchanged.

## 2026-11-16 — First public export
- Exported from the private development repository with `scripts/62_public_repo.py` (code, published papers' claims,
  schema, rating rules and docs; private material left out or redacted). Release: 0.9.1-preview.

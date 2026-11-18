---
name: claim-replication
description: Build or extend a claim-replication project for a field (a disease, tissue, model or question) - extract the published conclusions of a set of papers as structured claims, harmonise public datasets, re-test every claim at donor level in its own data and in every other dataset, rate it (verdict + certainty), and publish a claim store, knowledge graph, atlas and MCP server. Use when asked to "check which published findings hold up across datasets", "replicate the single-cell literature of X", "start a claim store for a new field", or to add papers, datasets or curated signatures to such a project.
---

# Claim replication: published claims -> re-tested, rated evidence

The method behind the reference project's claim stores (two fields so far), written so it can be run on another field.
It needs the `claimstore` engine and the project template of the reference repository (`template/`,
`scripts/63_new_project.py`); everything field-specific is config, one loader per dataset and one signature script per
paper. This skill lives in the reference repository's `.claude/skills/claim-replication/`; every project made from the
template gets a copy (update it from there, like the engine).

## What one project produces
- `claims/<paper>.yaml`: every main conclusion of each paper as a structured, testable claim with a short quote.
- Harmonised datasets: one schema, the same cell classes, donor-level pseudobulk, an edgeR DE cache per contrast.
- For every claim and dataset a re-test with an outcome (reproduced / direction_only / not_reproduced / contradicted
  / not_testable, with the reason), then a pooled estimate across datasets.
- A rating per claim: verdict (Supported, Partly supported, Inconclusive, No effect, Contradicted, Untested) and
  certainty (High, Moderate, Low, Very low) with the reasons (`claims/RATING.md`).
- A SQLite claim store, a knowledge graph, an atlas page, an MCP server, and a public release.

Calibration: the second field (8 papers, 204 claims, 2 datasets) took about 4 h of hands-on work and 3.5 h of
compute, with 439 lines of field code (two loaders, six signature scripts); everything else was config.

## When it applies
| Situation | Fit |
|---|---|
| Human tissue, single-cell / single-nucleus RNA, public datasets with donor metadata for the papers' contrasts | Direct: config, loaders, signature scripts |
| Other organism (mouse, organoids) | Direct with unit = animal or line; non-human evidence lowers certainty (indirectness) if the claims are about humans |
| Bulk RNA-seq or pseudobulk-only releases | Mostly: one class `ALL`; composition and cell-state claims become not testable |
| ATAC, spatial, proteomics, imaging | Claims, store and rating transfer; each needs a new test template in `claimstore.retest` and a harmonisation schema |
| No public data that can test the claims | Only the claims, store and atlas; every claim stays Untested |

## Decisions that belong to the user (ask before acting)
1. The field and the question (what counts as a case, a control, a stage).
2. The papers (inclusion rule and list) and the datasets.
3. Every download: list exact files, source and size first (`cli.py supplements --list`, the data portal's listing).
4. Case / control definitions and donor filters per dataset (e.g. "AD = ADNC Intermediate or High"; age >= 60).
5. What stays private (unpublished data or papers), the licence, and where it is published (repository, Pages, Zenodo).
6. Every commit, push, release and public page.
Record each decision with its date (project config comments, `EFFORT.md` "Analysis choices to review").

## Steps

### 1. Scope and project
- Make the project: `python3 -I scripts/63_new_project.py DIR --name ID --title TITLE --package PKG [--git]` in the
  reference repository. Its checks must pass at once.
- Fill the field's facts: `config/schema.toml` (region and cell-type vocabularies; keep the other vocabularies so
  claims stay comparable across fields), `config/contrasts.toml` (case / reference groups as obs filters; `[stage]`
  for a continuous stage, with a contrast code and label), `config/project.toml` (`[project]`, `[claims]`,
  `[tests]`, `[pipeline]`, `[release]`, `[atlas]`).
- Start `EFFORT.md`: one row per step with wall time, lines of code and what changed.

### 2. Papers and datasets
- Papers: the field's main single-cell papers with full text (PMC); rows in `inventory/papers.tsv` with `pmid`, `doi`,
  regions, modality, cohort and `reuses_data_of` (papers re-analysing another paper's data are not independent).
- Datasets: public, with raw counts and donor-level metadata for the claims' contrasts (diagnosis, staging, age, sex,
  PMI). Rows in `inventory/datasets.tsv`. Prefer datasets that cover the claims' regions and cell types; note brain
  bank or cohort, since datasets from one bank may share donors (reported, not merged).
- Download only after approval; check sizes against the published ones.

### 3. Claims (one file per paper)
Follow `references/extraction.md`. In short: every headline and main conclusion becomes one claim with the fields of
`claims/SCHEMA.md`; a verbatim quote of at most 15 words with its location; the comparison in the authors' words and
its orientation; the direction; the genes the claim names and where its full gene list is (`gene_source`); the test
template. Agents can extract in parallel (two papers per agent, full text from PMC); then merge and validate:
`cli.py merge` (claims.yaml, validation.tsv), `pytest tests`.

### 4. Contexts: which test, on which contrast, in which cells and regions
- `config/contexts.toml` maps each claim's comparison text to a contrast code (ordered regex rules), claim cell types
  to cell classes, claim regions to dataset regions, and flags no-difference claims (`null_pattern`). `cli.py contexts`
  writes `claims/contexts.tsv`; check that no claim is left without a contrast.
- Review the orientation of every claim against its contrast before testing (see Rules). Fix mismatches per claim in
  `claims/contrast_overrides.tsv` (claim_id, contrast, note); use `not_available` when no contrast in the data matches
  (an alternative trajectory, a measure the data lack).

### 5. Data: loaders, harmonisation, pseudobulk, DE
- One loader per dataset, `src/PKG/build/<dataset>.py`, `load()` returning an AnnData with raw integer counts and the
  obs columns of `config/harmonise.toml` (donor, sample, region, disease, sex, age, pmi, stage, author_cell_type),
  plus a map from the authors' labels to the project's cell classes (`AUTHOR_TO_CLASS`). Select donors here, with the
  rule in the loader's docstring. Large h5ad: `claimstore.h5io.read_csr` / `read_csr_many` (rows streamed into one
  matrix, no second copy).
- `cli.py harmonise DATASET` (loader, schema check, cell classes from the authors' labels checked against marker
  genes, donor table), `cli.py pseudobulk DATASET` (counts per donor x region x class and the unit table),
  `cli.py de-cache DATASET` (edgeR per region x class x contrast, and along the stage).
- Checks: nuclei and donors per group as expected; label agreement with markers (0.96-0.98 in the second field); every contrast of
  `contrasts.toml` fitted where the metadata allow (`de_log.tsv`); a run's memory peak noted for the next dataset.

### 6. Signatures (gene lists for gene-based claims)
- `cli.py supplements --list` lists each paper's supplementary files with sizes (nothing downloaded); after approval
  `cli.py supplements --only REGEX PAPER...`. Results hosted elsewhere (a code repository, a portal) go to
  `raw/supplements/<paper>/` with a `SOURCE.txt` naming the address and the commit or version.
- One script per paper reads its tables and writes `claims/signatures/<paper>.tsv` (rules in `claims/SIGNATURES.md`):
  DE lists split into de_up / de_down with the authors' threshold; state markers top 100; module or pathway members
  when listed. Set `[tests] signatures_per_class = true` and give each row its cell class, so a neuron list is not
  tested in glia.
- Validate each table against a number the paper reports (count of DEGs, genes shared by cell types) before using it.
- Claims whose test would not test the claim go to `<paper>_unresolved.tsv` with a reason (see Rules).
- `cli.py signatures` merges curated tables with the genes the claims name; a curated signature replaces the named
  genes of the claims it covers.

### 7. Re-tests and meta-analysis
- For each dataset: `cli.py tests internal DATASET` (claims whose source data it is) and `cli.py tests cross DATASET`
  (every other claim, regions retargeted to the dataset's; claims whose own data are not in the project included),
  then `cli.py meta`.
- Read the outcome counts and the reasons for not_testable (`note`): they say what to fix next (missing signatures,
  a contrast the data lack, too few donors, a template the engine does not have).

### 8. Store, rating, atlas, server
- `cli.py store`, `cli.py rating`, `cli.py store` again (loads the ratings), `cli.py atlas`, `cli.py serve`.
- After every change, diff the ratings against the previous version (git) and explain every verdict that moved.

### 9. Release
- Public release: `cli.py release VERSION`; public atlas `cli.py atlas --public --db <release store>`; audit that no
  private dataset, paper or name appears; public repository with a fresh history if the working history holds private
  material; Pages from `docs/`; Zenodo for the store, under the user's account.

### Optional layers (the reference project has them)
Findings of the project's own analyses rated like claims (`claims/findings.yaml`), consensus statements over groups of
claims, reference knowledge from other evidence types (`claims/background.yaml`), and the question loop (atlas
questions -> new analysis -> finding; the `atlas-questions` skill).

## Rules (hold in every field)
- **Donor-level units.** Every test compares donors (pseudobulk per donor x region x class), never cells.
- **Test everywhere.** Each claim in its own data (internal replication) and in every other dataset (cross).
- **Untested is not false.** It means no dataset could test it; report it with the reason.
- **Never force a test.** If a re-test of the genes would not test what the claim says (overlap or specificity of gene
  sets, ranking between subpopulations, order along a trajectory, a pathway named only in a figure, a comparison
  between two studies), record the claim as unresolved with the reason instead.
- **Orientation.** The direction is relative to the contrast's case group. A positive association with *decline* is
  negative along *better cognition*; a trajectory with no counterpart in the data is `not_available`, not the nearest
  contrast. Check each claim whose comparison mentions a second variable.
- **No-difference claims are claims.** They are rated on equivalence (`Supported (no effect)`); make sure the
  `null_pattern` catches the field's phrasing ("do not differ", "does not change", "no ... changes").
- **Scope.** A claim holds for the regions, cell types and contrast tested; other-region and proxy evidence lowers
  certainty (indirectness), evidence from the source cohort counts as dependent.
- **Uniform curation.** The same selection rule for every paper; record the rule in every signature row.
- **Quotes of at most 15 words**, each with its location. No reproduction of tables or text beyond that.
- **Private stays private.** Unpublished data and papers are listed in `[release]` and filtered from every export.
- **Respect refusals.** If a server refuses scripted downloads (403, a bot check), stop and ask; use another
  legitimate source (the publisher) only with approval.
- **Engine changes keep the other fields identical.** Prove it (old vs new outputs byte-identical) or make the change
  opt-in through config.

## Pitfalls seen in the first two fields
- A no-difference pattern written for one field missed every no-difference claim of the next.
- Comparison text mentioning cognition sent disease-progression claims to a "better cognition" contrast (wrong sign).
- Claims tested on 3-5 named genes changed verdict when the paper's full table was used (Contradicted -> Supported
  and Supported -> Inconclusive both happened): curate signatures before trusting gene-level verdicts.
- Per-cell-type DE tables applied as one gene list tested neuron genes in glia.
- Region codes from a brief were wrong (Brodmann area 10 is frontopolar, not DLPFC): keep a glossary.
- Supplementary files of author manuscripts are refused by PMC and Europe PMC to scripts; tables on external buckets
  are only links in the supplement.
- Large datasets (> 2^31 stored values) need 64-bit indices and in-place reading, or the loader runs out of memory.
- Release metadata can lack the staging a contrast needs (one release had no neuropathological stage): those tests are not testable.

## Effort log and write-up
Keep `EFFORT.md`: per step the wall time, compute time, lines of code, config changed, and what the field forced back
into the engine or template; the results table (verdicts before and after each change); why claims stay untested
(counts by reason); the analysis choices to review. It is the record a methods section and a reviewer need.

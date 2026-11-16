# Claim store — design proposal (draft, 2026-11-13)

Goal: one store where every claim (published, from our re-analyses, or raised by a question) carries its full context and provenance, is linked to the datasets it comes from and was tested in, and is linked to other claims, so that it can be browsed, queried by an AI, and extended.

## 1. What we have now, and what is missing

| piece | where | problem |
|---|---|---|
| claims (456) | `claims/<paper>.yaml` → `claims.tsv` | one flat row; context fields are free text |
| papers, cohorts | `inventory/papers.tsv`, `datasets.tsv` | cohort facts (n, regions, access) are free text; not joined to claims by key |
| contrasts | `claims/contexts.tsv`, `contrast_overrides.tsv` | normalised contrast codes exist only for testing; `comparison` in claims has 155 distinct strings |
| test results | lustre `results/tests/{internal,cross}/*.tsv`, `results/meta/*` | the claim × dataset links exist, but outside git and only as engine output |
| consensus | `consensus_statements.tsv` | claim links as `;`-joined strings; no relation types |
| our own findings | manuscript text, `CHANGELOG.md`, `results/subgroups/` | not stored as claims at all (axis 2, clock, specificity, ...) |

Missing context, concretely:
- **Organism / tissue / sample type** are implicit. "Human post-mortem brain" is a project-level assumption, so a claim read alone ("microglia expand in PD") loses it. Some claims are not human at all (kamath2022-C02 is cross-species, C03 is macaque Slide-seq).
- **Cohort**: ancestry is stored only for 3 datasets (lee2026 `genetic_ancestry`, prashant2024 `ancestry`, wang2024 `race`); brain bank, age, sex balance, disease definition and Lewy stage of the donors are not attached to claims.
- **Region**: 11 coarse values, 62 claims are `multiple`, `striatum` vs `putamen` mixed, sub-region (SNpc ventral tier, cortical layer, grey vs white matter) only in free text; no ontology ID.
- **Comparison / disease definition**: "PD/LBD vs control" hides who is in each group (PDD? iLBD? Braak range?).
- **Sampling**: FANS-enriched nuclei (kamath2022 NR4A2 sort), grey-matter-only dissections, etc., change what a composition claim can mean, but are in notes only.
- **Scope vs evidence**: the store does not separate *where the authors observed it* from *where it has been tested since*.

## 2. Model: a claim is an assertion + context + provenance + evidence

The model is close to a nanopublication or a Biolink association, kept relational and simple.

```
            ┌──────────── Context (reusable, ontology-coded) ────────────┐
            │ organism  sample_type  region(s)  cell type/state  assay  │
            │ cohort    disease/contrast definition  sampling notes     │
            └──────────────────────────────┬────────────────────────────┘
                                           │ scoped by
 Source ──asserts──▶ CLAIM  ── subject ─ predicate ─ object ── direction
 (paper / our analysis /       │                                  
  user question)               ├── derived_from ──▶ DATASET (role: source)
                               ├── tested_in ─────▶ DATASET (outcome, effect, n, region match)
                               └── relation ──────▶ CLAIM (supports / contradicts / refines /
                                                            same_as / member_of consensus)
```

### Core entities (tables)

| table | key fields |
|---|---|
| `source` | id, kind (`paper` / `our_analysis` / `question`), citation or script + commit, date, author |
| `cohort` | id, brain bank(s), n donors by group, **ancestry composition** (or `not_reported`), age range, sex ratio, PMI, disease definitions, Lewy staging scheme |
| `dataset` | id, cohort, assay (EFO), regions (UBERON), enrichment/sorting, n nuclei, accession, access, harmonised path |
| `donor_group` | dataset, group label, definition (diagnosis, Braak range), n — so a contrast points at real groups |
| `entity` | cell types (CL), cell states (our ids + author names + signature id), genes (HGNC), pathways (GO/Reactome), programmes (axis 2, heat-shock) |
| `claim` | id, source, statement (verbatim-ish), **subject** (entity), **predicate** (controlled), **object/contrast** (group A vs B), direction, claim_type, centrality, unit of analysis, reported stat, context id |
| `context` | organism (NCBITaxon:9606), sample type (post-mortem brain), region(s) + sub-region, cell type/state, assay, cohort/dataset, contrast — one row per distinct scope, reused by many claims |
| `claim_dataset` | claim, dataset, role (`source` / `tested`), outcome, effect, CI, p, n case/ref, region match, test template, run id |
| `claim_relation` | claim A, relation, claim B, basis (`curated` / `computed`), note |
| `question` | id, text, asked by, date, status; links to the analyses and claims it produced |

Controlled predicates (small set, mapped from `claim_type` + `direction`): `increased_in`, `decreased_in`, `exists_in`, `enriched_in`, `depleted_in`, `correlates_with`, `heritability_enriched_in`, `accessible_in`, `colocalises_with`, `not_different_in` (null claims, already treated by the equivalence rule).

### Two levels of statement

- **Short statement**: as now, for display.
- **Contextualised statement**, generated from the fields, never typed by hand, so it cannot drift. Example for kamath2022-C06:

> In **human post-mortem** substantia nigra pars compacta (**NR4A2-sorted nuclei, snRNA-seq**), the SOX6_AGTR1 dopaminergic subtype is **depleted** in PD/LBD (10 donors) versus controls (8 donors) [cohort: Kamath 2022, ancestry not reported; unit: cells, MASC]. *Since tested in 6 datasets: 1 reproduced, 5 direction only, 0 contradicted.*

That sentence is what an AI reads or embeds, and what the browser shows on hover.

### Context inheritance

Most context lives on the dataset and cohort, not on the claim. A claim stores only what *narrows* it (sub-region, cell state, a specific donor group). The full context is resolved claim → context → dataset → cohort at build time. Ancestry is then entered once per cohort (about 20 rows) instead of 456 times, and it is explicit when it is unknown.

## 3. Storage: three layers

1. **Source of truth: YAML/TSV in git** (as now, extended). It is human-editable, reviewable in diffs and agent-friendly. `claims/<source>.yaml` gains `subject`, `predicate`, `contrast` (structured), `context` overrides; new `inventory/cohorts.yaml` with ancestry etc.; `claims/relations.tsv`; `questions/<id>.yaml`.
2. **Compiled store: one DuckDB (or SQLite) file**, built by a script from layer 1 plus the lustre test outputs. It has the tables above, with validation (every region a known UBERON id, every claim with a source dataset, every contrast pointing at donor groups). This is what everything downstream queries.
3. **Views**, all generated from the compiled store:
   - **Browser**: two-way pages, claim → source dataset + every test (outcome, effect, region match), and dataset → claims derived from it + claims tested in it, with filters (region, cell type, ancestry, outcome). Options: extend the existing `claim_browser.html` (static, shareable as an Artifact), or Datasette over the SQLite (zero code, local faceted browsing and SQL).
   - **AI access**: (a) a JSONL "knowledge pack" of contextualised statements + evidence for retrieval, and/or (b) a small MCP server exposing `search_claims`, `get_claim`, `claims_for_dataset`, `datasets_for_claim`, `add_question` over the DB, so Claude can answer from the store and cite claim ids instead of from memory.
   - **Graph export** (optional, later): the same tables map to nodes and edges (claims, entities, datasets) for Neo4j / Cytoscape if network views are wanted. At ~500 to a few thousand claims a graph database is not needed for storage.

Why not a graph DB or a single big JSON as the primary store: diffs and review get hard, and the relational form already supports every query we need (claim ↔ dataset ↔ cohort ↔ entity). Why not only TSV: no referential checks, and joins across 6 tables become ad hoc pandas code in every script.

## 4. Claims from our analyses and from questions

Same table, `source.kind = our_analysis | question`, so the paper's findings (axis 2 is whole-brain, astrocyte age gap +3.5 y, CIRBP loss shared with AD/FTD, ...) and side-project results become first-class claims. Each gets:
- provenance: script, commit, results file, date;
- `claim_dataset` rows for every dataset it used;

Question loop: question → `question` row → analysis (script) → new claims with `source = question:<id>` → relations to existing claims. The store then also records *why* each analysis was run.

Status field for every claim: `asserted` (as published) → `tested` → `consensus member`; plus `superseded_by` for our own claims when an analysis is rerun.

## 5. Ontologies (IDs stored next to labels, labels stay human)

organism NCBITaxon; region UBERON (e.g. SNpc UBERON:0001965, putamen UBERON:0001874, DMNX UBERON:0002870); cell type CL; disease MONDO (PD MONDO:0005180, DLB MONDO:0007488); assay EFO; ancestry HANCESTRO; genes HGNC. (IDs above are written from memory; check them against OLS when the vocabularies are built.) This also makes it possible to match to CELLxGENE and other external stores later.

## 6. Proposed order of work

1. `inventory/cohorts.yaml`: one entry per cohort with ancestry, brain bank, age, sex, disease definitions, staging (from papers, supplements and harmonised obs; `not_reported` stays explicit).
2. Region/cell type/disease vocabularies with ontology IDs; map the existing values (mechanical for most, about 60 `multiple` claims by hand).
3. Structured `contrast` + `subject/predicate` for existing claims: derive from `contexts.tsv` and `contrast_overrides.tsv` where possible, then review the rest.
4. Build script → DuckDB with validation; import lustre test outputs as `claim_dataset`.
5. Contextualised-statement generator; two-way browser (dataset ↔ claim).
6. Our own analysis claims (manuscript Results + side projects) entered as `our_analysis`.
7. AI access: knowledge pack, then MCP server; `question` workflow.

Steps 1–5 are a working prototype; 6–7 make it the place where new questions land.

## 7. Prototype status (2026-11-13)

Built: steps 1–5 (`scripts/51a_vocab.py`, `51b_claim_store.py`, `51c_store_browser.py`; outputs in `claims/vocab/`, `claims/store/`, `results/store/claim_store.sqlite`, `reports/pd_atlas.html`). The browser view is a visual atlas: a landing summary of what changes in PD (charts from results/biology) with drill-down pages per cell type, programme, region, consensus statement, claim, dataset and gene. Differences from the plan above:
- SQLite, not DuckDB (not installed in the env; SQLite also works directly with Datasette).
- Structure is derived, not hand-entered: `claim_structure.tsv` is generated from the existing claim fields and `contexts.tsv`; hand fixes go in `claims/structure_overrides.tsv`. 62 claims with region `multiple` and 15 without a contrast are flagged `needs_review`.
- Cohort context comes from an extracted draft until `inventory/cohorts.yaml` is reviewed.
- Ontology IDs in section 5 were checked against OLS4: the snRNA-seq ID is EFO:0009809 and substantia nigra is UBERON:0002038 (SNpc UBERON:0001965 is not used, as datasets dissect whole SN).
Not built yet: step 6 (our own analysis claims) and step 7 (MCP server, question workflow).

Update (later 2026-11-13): step 6 and the question workflow are built (`claims/findings.yaml`, `claims/questions.yaml`, scripts 53-56). Step 7 is now built too:
- **Knowledge graph** (`src/pdrep/kg.py`, export `scripts/57_knowledge_graph.py`): the store read as typed nodes (claim, paper, dataset, cell type, region, gene, theme, consensus statement, reference knowledge, question) and relations (asserts, about_cell_type, in_region, names_gene, has_theme, derived_from, tested_in with outcome, supports/opposes_statement, agrees_with, establishes, bears_on, supports, consistent_with, contrasts_with, concerns_region, generated, covers_region, cites, produced). About 1,130 nodes and 8,350 edges; `not_testable` re-tests and signature member genes are left out. Exports: `claims/store/kg_nodes.tsv`, `kg_edges.tsv` (git), `results/store/knowledge_graph.graphml` (Cytoscape, Gephi) and `.json` (networkx node-link). The relational store stays the source of truth; the graph is rebuilt from it.
- **Atlas Graph page** (`#graph`): an overview map (consensus statements, reference knowledge and the cell types their claims are about, weighted by claim count) and a neighbourhood view for any record (`#graph/<node id>`, 1 or 2 steps, optional re-tests, genes, themes); claim and statement pages link to it.
- **3D view** (`#graph3d`): claims, consensus statements and reference knowledge full screen, navigable in 3D, laid out by shared context (papers, datasets, genes, cell types, regions and themes act as invisible anchors) or by meaning (embedding neighbours, `kg.similarity_links`); every other category is a colouring.
- **MCP server** (`src/pdrep/mcp_server.py`, registered in `.mcp.json`): stdio JSON-RPC, standard library only. Tools: overview, search, search_claims, get_claim, get_statement, get_dataset, gene, neighbours, connect, schema, sql (read-only), list_questions, add_question (the only write: appends to `claims/questions.yaml`; disabled with PDREP_MCP_READONLY=1).
- **Remote server** (`src/pdrep/mcp_http.py`, deployment in `deploy/mcp/`): the same tools for any AI system, as MCP over Streamable HTTP (`/mcp`) and as REST + OpenAPI (`/api/tools/<tool>`, `/openapi.json`). Read-only (add_question off), SQL time limit, row cap, per-IP rate limit, optional bearer token; container image with a snapshot of the store.

# Claim signatures

Turns each claim that names a cell state, genes or a pathway into the gene list the re-test needs. One script per paper,
`scripts/signatures/<paper>.py`, reads that paper's supplementary files from
`/lustre/groups/ml01/datasets/projects/2025_pdstriatum_sergio.salas/pd_sc_replication/raw/supplements/<paper>/`
and writes `claims/signatures/<paper>.tsv`. Run with `~/miniconda3/envs/rapids_singlecell/bin/python -I`.

## Output columns (`claims/signatures/<paper>.tsv`, long format, one row per gene)

| column | meaning |
|---|---|
| `signature_id` | `<paper>__<short_name>`, e.g. `kamath2022__SOX6_AGTR1_markers`, `wang2024__astro_DEG_up` |
| `paper` | paper id |
| `claim_ids` | claim ids that use this signature, `;`-separated |
| `kind` | `state_markers` (genes defining a subpopulation) / `de_up` / `de_down` (PD-vs-control or other comparison) / `pathway_genes` / `named_genes` (genes named in the claim text only) |
| `cell_type` | broad class (schema vocabulary) |
| `cell_state` | authors' state name, or empty |
| `comparison` | for de_*: the comparison; empty for markers |
| `gene` | HGNC symbol as given by the authors (do not map or rename) |
| `rank` | 1 = strongest, by the authors' ordering statistic |
| `stat` | value of that statistic (logFC, score, ...) |
| `padj` | adjusted p as reported, empty if not reported |
| `source` | `file|sheet` or `text` |
| `selection` | the rule used, e.g. `padj<0.05 & logFC>0.25, top 100 by logFC` |

## Selection rules (apply uniformly)

- **State markers:** genes the authors list as markers for that state; keep those with reported adj. p < 0.05 and positive logFC when the stats exist, ranked by logFC (or the authors' score); cap at the top 100. If only a short named list exists, use it and set `selection: text list`.
- **DE gene lists:** all genes the authors call significant for that comparison and cell type, split into `de_up` / `de_down`, with their threshold in `selection`. No cap (the full table is used for direction concordance).
- **Named genes:** when a claim rests on genes named in the text only (e.g. "HSPA1A is up in astrocytes"), emit one `named_genes` signature with those genes and `source: text`.
- **Pathways:** if the authors give the member genes (leading edge or table), use them; if only a pathway name is given (GO/KEGG/Hallmark id), emit no genes and record the term in `claims/signatures/<paper>_unresolved.tsv` (claim_id, term, reason) so it can be resolved from MSigDB centrally.

Every claim with `test_template` in state_signature, geneset_shift, de_pseudobulk, or composition with a non-empty `cell_state` must end up either in a signature (`claim_ids`) or in `<paper>_unresolved.tsv` with a reason (e.g. "marker table not in supplement").

# Claim signatures

Turns each claim that names a cell state, genes or a pathway into the gene list its re-test needs.

1. `python3 -I src/claimstore/cli.py supplements --list [PAPER...]` lists each paper's supplementary files (from its
   PMC full-text XML) with sizes where the server reports them; nothing is downloaded. Get the downloads approved,
   then run it without `--list`, with `--only REGEX` to keep just the approved files. Files go to
   `raw/supplements/<paper>/`; `inventory/supplements_manifest.tsv` records what came from where. Results the authors
   host elsewhere (a code repository, a data portal) go to `raw/supplements/<paper>/` too, with a `SOURCE.txt` naming
   the address and version.
2. One script per paper (for example `src/{{package}}/signatures/<paper>.py`) reads those files and writes
   `claims/signatures/<paper>.tsv`, plus `claims/signatures/<paper>_unresolved.tsv` for what it could not resolve.
3. `python3 -I src/claimstore/cli.py signatures` merges them into `claims/signatures_all.tsv.gz`; a curated signature
   replaces the genes its claims name in their text. Then re-run the re-tests (`tests`, `meta`, `store`, `rating`).

## Output columns (`claims/signatures/<paper>.tsv`, long format, one row per gene)

| column | meaning |
|---|---|
| `signature_id` | `<paper>__<short_name>`, e.g. `<paper>__microglia_state3_markers`, `<paper>__astro_DEG_up` |
| `paper` | paper id |
| `claim_ids` | claim ids that use this signature, `;`-separated |
| `kind` | `state_markers` (genes defining a subpopulation) / `de_up` / `de_down` (case vs reference, or the claim's own comparison) / `pathway_genes` / `named_genes` (genes named in the claim text only) |
| `cell_type` | broad class (the schema's cell-type vocabulary) |
| `cell_state` | the authors' state name, or empty |
| `comparison` | for de_*: the comparison; empty for markers |
| `gene` | HGNC symbol as given by the authors (do not map or rename; the re-tests map aliases) |
| `rank` | 1 = strongest, by the authors' ordering statistic |
| `stat` | value of that statistic (logFC, score, ...) |
| `padj` | adjusted p as reported, empty if not reported |
| `source` | `file|sheet` or `text` |
| `selection` | the rule used, e.g. `padj<0.05 & logFC>0.25, top 100 by logFC` |

## Selection rules (apply uniformly)

- **State markers:** genes the authors list as markers for that state; keep those with reported adjusted p < 0.05 and
  positive logFC when the statistics exist, ranked by logFC (or the authors' score); cap at the top 100. If only a
  short named list exists, use it and set `selection: text list`.
- **DE gene lists:** all genes the authors call significant for that comparison and cell type, split into `de_up` /
  `de_down`, with their threshold in `selection`. No cap (the full table is used for direction concordance).
- **Named genes:** when a claim rests on genes named in the text only, the `signatures` command already uses the
  claim's `genes` field; a script adds nothing.
- **Pathways:** if the authors give the member genes (leading edge or table), use them; if only a pathway name is
  given (GO/KEGG/Hallmark id), emit no genes and record the term in `<paper>_unresolved.tsv` (claim_id, term, reason).

Every claim with `test_template` in state_signature, geneset_shift, de_pseudobulk, or composition with a non-empty
`cell_state` should end up either in a signature (`claim_ids`) or in `<paper>_unresolved.tsv` with a reason (e.g.
"marker table not in supplement").

# Claim schema

One file per paper: `claims/<paper_id>.yaml`, a YAML list of claims. Merged into `claims/claims.yaml` by `scripts/05_merge_claims.py`.

## What counts as a claim

A main conclusion the authors draw about PD (or about the PD brain versus controls) from their single-cell / spatial data: what the abstract, the Results section headings and the figure titles assert. Typically 5–15 per paper. Not counted: methods descriptions, atlas statistics ("we profiled 200,000 nuclei"), conclusions drawn only from mouse, iPSC or bulk data (record these only when the paper presents them as validation of a single-cell claim, in `notes`).

Split compound sentences into separate claims when they can be true or false independently ("astrocytes and microglia are depleted" = 2 claims). Keep the authors' cell-type and state names.

## Fields

| field | values / meaning |
|---|---|
| `id` | `<paper_id>-C<nn>` (C01, C02, ... in order of appearance) |
| `paper` | paper id |
| `statement` | the conclusion in one sentence, close to the authors' wording |
| `quote` | short verbatim quote supporting it (≤ 40 words) |
| `location` | figure / panel / section, e.g. `Fig. 3c; Results "Loss of SOX6+ ..."` |
| `centrality` | `headline` (in abstract or title) / `main` (Results heading or main figure) / `secondary` (supplementary or passing) |
| `modality` | snRNA / scRNA / multiome-RNA / snATAC / multiome-ATAC / spatial-Visium / spatial-Xenium / spatial-GeoMx / spatial-other / imaging |
| `region` | region(s) the claim is about (use: SN, VTA, midbrain, putamen, caudate, striatum, PFC, DLPFC, cingulate, temporal, MTG, motor, visual, parietal, amygdala, hippocampus, cerebellum, DMNX, GPi, multiple) |
| `cell_type` | broad class: DA neurons, Neurons (other), Excitatory neurons, Inhibitory neurons, MSNs, Astrocytes, Microglia, Oligodendrocytes, OPCs, Endothelial, Pericytes, Fibroblasts, T cells, Immune (other), All |
| `cell_state` | named subpopulation or signature (e.g. `SOX6_AGTR1`, `GPNMB-high microglia`, `UPR program`), or empty |
| `claim_type` | composition / de_gene / de_pathway / state_existence / state_enrichment / genetic_enrichment / communication / pathology_association / trajectory / chromatin / spatial_colocalisation / clinical_association / other |
| `direction` | up / down / enriched / depleted / exists / associated_pos / associated_neg / NA |
| `comparison` | e.g. `PD vs control`, `PDD vs PD`, `Braak 5-6 vs 3-4`, `within PD` |
| `genes` | genes needed to test the claim (list), as reported in the text/figure |
| `gene_source` | where the full gene list is (e.g. `Supplementary Table 4, sheet "astro_DEG"`), or `text` |
| `reported_stat` | effect size / statistic / p or FDR as reported, and n |
| `unit_of_analysis` | cells / donors / samples / ROIs — what the authors' test treated as independent |
| `method` | test used by the authors (e.g. `Wilcoxon on cells`, `MAST + donor RE`, `scCODA`, `MAGMA`) |
| `testable_with` | rna / atac / spatial / metadata / not_testable |
| `test_template` | composition_donor / de_pseudobulk / geneset_shift / state_signature / genetic_enrichment / communication_liana / pathology_correlation / trajectory / chromatin_da / spatial_neighbourhood / not_testable |
| `requires` | metadata or features needed beyond PD/control labels (e.g. `Braak stage`, `disease duration`, `sorted DA neurons`) |
| `notes` | caveats, ambiguities |

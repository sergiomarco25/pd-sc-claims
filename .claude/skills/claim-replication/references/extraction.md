# Claim extraction brief

How to turn a paper into `claims/<paper_id>.yaml`. Fields and vocabularies: `config/schema.toml`; format:
`claims/SCHEMA.md`; example: `claims/example2024.yaml`. Give this brief, the schema and one finished paper file to each
extraction agent (two papers per agent works well).

## 1. Get the full text
- PMC full text: Europe PMC `https://www.ebi.ac.uk/europepmc/webservices/rest/<PMCID>/fullTextXML`. For NIH author
  manuscripts outside the open-access subset it fails; NCBI `efetch.fcgi?db=pmc&id=<number>` or the PMC OAI-PMH service
  return the same JATS XML. PMCIDs from PMIDs: the NCBI ID converter.
- Read the abstract, every Results section, figure legends and the Methods (groups, unit of analysis, statistics).
- If only the abstract is available (paywalled, not in PMC), extract from it, say so in `location` and `notes`, and do
  not invent figure panels or statistics.

## 2. What is a claim
- One conclusion the authors draw from their own data, specific enough to test in another dataset: a group
  difference, an association, a subpopulation and its markers, a gene's change in a cell type, a pathway's change, a
  relation to pathology or clinical measures. Results of re-analysed public data count, marked in `notes`.
- One claim per testable statement. Split "X is up in astrocytes and down in microglia" into two claims. Merge
  restatements of the same result (abstract, Results, Discussion) into one claim, citing the main location.
- Not claims: background from other studies, methods, speculation in the Discussion, technical validations without a
  biological conclusion.
- Cover every headline and main conclusion; add secondary ones when they are testable. Typical yield: 20-35 per paper.

## 3. Fields
| field | rule |
|---|---|
| `id` | `<paper_id>-C<nn>`, in order of appearance |
| `statement` | the conclusion in your own words, one sentence, self-contained (cell type, region, groups, direction) |
| `quote` | the authors' words supporting it, verbatim, at most 15 words |
| `location` | figure and panel, or section ("Fig. 3b; Results 'Microglial states'") |
| `centrality` | headline (abstract or title) / main (Results heading or main figure) / secondary |
| `modality`, `region`, `cell_type` | from the vocabularies; regions as the paper defines them (see 4) |
| `cell_state` | the authors' subpopulation or signature name, as written ("Mic.12", "a1; a6"), or empty |
| `claim_type` | composition, de_gene, de_pathway, state_existence, state_enrichment, pathology_association, clinical_association, trajectory, genetic_enrichment, chromatin, spatial_colocalisation, communication, other |
| `direction` | up / down / enriched / depleted / exists / associated_pos / associated_neg / NA; relative to the case group of `comparison` (see 5) |
| `comparison` | the groups or variable in the authors' terms: "AD vs control", "stage 6 vs stage 0", "along global cognition", "within AD donors" |
| `genes` | genes the claim names, as written |
| `gene_source` | where the full gene list is ("Supplementary Table 2, sheet DEGs", "authors' GitHub, gpath results"), or "text" |
| `reported_stat` | effect, statistic, p or FDR and n, as reported |
| `unit_of_analysis` | what the authors treated as independent: cells, donors, samples, ROIs |
| `method` | the authors' test |
| `testable_with` | rna / atac / spatial / metadata / not_testable |
| `test_template` | see 6 |
| `requires` | metadata the test needs beyond diagnosis (stage, APOE, sex, cognition scores) |
| `notes` | caveats: cell-level statistics, small n, re-analysed data, ambiguity, abstract only |

## 4. Regions
Use the project's region codes. Map the paper's description carefully and note the original wording in `notes`
when the mapping is not exact (Brodmann area 10 is frontopolar, not DLPFC; BA9/46 is DLPFC). When in doubt, ask and
record the decision in the project's region glossary.

## 5. Direction and orientation
- The direction is relative to the case group of the comparison as the authors state it ("higher in AD" -> up for
  "AD vs control").
- Associations: associated_pos / associated_neg with the variable as written. Write the variable so its sign is
  unambiguous: "along cognitive decline" and "along global cognition" have opposite signs. Say which in `comparison`.
- Trajectories, pseudotime and modelled programmes: name the trajectory or score ("along the prAD trajectory",
  "along CPS"); the contexts step decides whether any dataset has a matching contrast.
- **No-difference claims** ("do not differ", "unchanged", "no subpopulation changes"): direction NA, and write the
  absence explicitly in the statement so the null pattern catches it; mention it in `notes`.

## 6. Test template by claim type
| claim | template | needs |
|---|---|---|
| a cell type's share changes | composition_donor | cell classes in the data |
| a subpopulation's share changes | composition_donor with `cell_state` | the state's markers (gene_source) |
| gene(s) change in a cell type | de_pseudobulk | the genes |
| a pathway / gene set changes | geneset_shift | member genes (gene_source); a pathway name alone is not enough |
| a subpopulation exists / has markers | state_signature | markers |
| a state or score tracks pathology or clinical measures | pathology_correlation | the stage or measure in the data |
| genetic risk enriched in a cell type | genetic_enrichment | the field's genetic test |
| trajectory, chromatin, spatial, communication | trajectory / chromatin_da / spatial_neighbourhood / communication_liana | matching data |
| nothing in expression data could test it | not_testable | reason in `notes` |

## 7. Per paper, besides the claims
- A row in `inventory/papers.tsv` (pmid, doi, regions, modality, n cases / controls, nuclei, accession, data access,
  author labels available, cohort, `reuses_data_of`).
- A short note in `claims/notes/<paper_id>.md`: what was read, what was left out and why, where the supplementary
  tables are and what they hold (this feeds the signature step).

## 8. Checks before merging
- Every quote <= 15 words and found verbatim in the text; every location exists.
- No claim without comparison, direction and template; no-difference claims phrased explicitly.
- `cli.py merge` passes; `validation.tsv` lists no errors; ids unique and consecutive.

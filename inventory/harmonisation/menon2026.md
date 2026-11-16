# menon2026 (Banner Sun Health multiome atlas) — harmonisation notes

- Loader: `src/pdrep/build/menon2026.py`.
- Builds:
  - 41861183 was cancelled because of the `_io` overflow bug.
  - **41861545** is the final build: 500 GB requested, 117 GB peak, 38 min.
- Output: 4,502,568 nuclei × 25,183 genes, 101 donors, 505 donor × region units.

## Source
- `$D/external_validation/menon2026/menon2026_complete_rna.h5ad` (CELLxGENE collection 16876983-…, RNA modality).
- `X` is normalised. **Raw integer counts are in `raw/X`**, with the same 25,183 genes. These are streamed in row chunks; the logic follows revision `33_menon_counts_cache.py`.
- Symbols come from `feature_name`; Ensembl IDs are in `var.gene_ids`.
- Donor metadata comes from Suppl. Table 1, via `revision/results/30_menon2026/donor_metadata.csv` (read-only). All 101 BannerIDs match `donor_id`.

## obs
- `donor` = donor_id; `sample` = donor_region; `batch` = pool prefix of the barcode (e.g. X056_G8).
- `region`: SN, PUT → putamen, CING → cingulate, MTG, CBL → cerebellum.
- `disease`: Control = yes → control. PD = yes and `dementia_yn` = 1 → **PDD** (46). Other PD → PD (34).
  - This follows the HARMONISATION rule (use PDD when the source says dementia).
  - **The paper treats all 80 as PD**, so PD+PDD must be pooled to reproduce its contrasts.
  - `disease_raw` is the CELLxGENE label. `clinical_dx` (neurological_dx) and `dementia` are kept.
- `braak_lb` = NaN. The cohort is staged with the Unified Staging System for LB disorders, which is not Braak. This stage is kept as:
  - `unified_lb_stage`: 0, 2 = IIa, 3 = III, 4 = IV
  - `unified_lb_stage_raw`: the original string
- `braak_nft` comes from 'Braak score' (I–IV → 1–4).
- Other columns: `age` = age_at_death; `pmi` = PMI; sex from CELLxGENE; `cell_type_cellxgene`; `brain_bank`.

## Counts vs paper (claims/notes/menon2026.md)
| | paper | harmonised |
|---|---|---|
| donors | 101 (80 PD + 21 control) | 101: 34 PD + 46 PDD + 21 control |
| samples | 505 (5 regions/donor) | 505 units; every donor has all 5 regions |
| nuclei (good snRNA) | 4.5 M | 4,502,568 |
| cell types | 30 subtypes, 7 broad classes | 30 author labels |
| Unified LB stage | I 21, IIa 6, III 37, IV 36 | 0 ("No Lewy bodies") 21, IIa 6, III 37, IV 37 |

- The paper's "Stage I N = 21" is the 21 controls. Suppl. Table 1 gives them as "0. No Lewy bodies", which resolves the ambiguity in the notes.
- The paper's text gives IV = 36, while the table has 37 (sum 101).

Cells per class:

| class | nuclei |
|---|---|
| Oligodendrocytes | 1,656,143 |
| Excitatory neurons | 1,501,803 (incl. 849,359 cerebellar granule cells) |
| Astrocytes | 438,271 |
| Inhibitory neurons | 257,754 |
| Microglia | 247,942 |
| MSNs | 203,511 |
| OPCs | 100,041 |
| DA neurons | 43,647 |
| Endothelial | 19,838 |
| Immune (other) | 13,539 |
| Fibroblasts | 8,178 |
| Pericytes | 6,566 |
| T cells | 3,357 |
| Ependymal | 1,978 |

## Label mapping
`author_cell_type` = obs.author_cell_type (30 labels).

| author labels | class |
|---|---|
| Oligo | Oligodendrocytes |
| OPC | OPCs |
| Astrocyte, Bergmann-Glia | Astrocytes |
| Microglia | Microglia |
| PV-Macrophage | Immune (other) |
| T-cell | T cells |
| Excitatory, RORB-Neuron, Granule (cerebellar, glutamatergic) | Excitatory neurons |
| DRD1/DRD2/DRD1-DRD2-MSNs | MSNs |
| Dopaminergic | DA neurons |
| VIP, PVALB, SST, LAMP5, Chandelier, GABA-Put, GABA-SN, MLI, Golgi, Purkinje | Inhibitory neurons |
| Capillary | Endothelial |
| Pericyte, SMC | Pericytes |
| Fibroblast, PV-Fibroblast | Fibroblasts |
| Ependymal | Ependymal |

## Label check (original harmonise.label_check)
Cell-weighted agreement is 0.548.
- **All neuron labels are below 0.5**, with marker majority "Neurons (other)". This includes Dopaminergic (0.05), the MSN labels (0.18–0.27), Excitatory (0.11), Granule (0.001) and most interneuron labels. This is the pan-neuronal-score artefact; re-run with the revised check (`10_relabel.py`).
- **Other labels below 0.5:**
  - MLI: 0.39, majority OPCs.
  - Ependymal: 81 cells checked, 0.0, majority Endothelial.
  - Fibroblast: 139 cells checked, 0.0, majority Oligodendrocytes.

  These are small populations; the author labels are kept.
- Glia and vascular labels are 0.52–0.98.

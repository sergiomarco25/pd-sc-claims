# prashant2024 (AMP-PD multiregion atlas) — harmonisation notes

- Loader: `src/pdrep/build/prashant2024.py`.
- Builds:
  - 41861083 failed: int32 overflow in the loader's CSR stacking, since fixed in `build/_io.py`.
  - **41861544** succeeded: 450 GB requested, 131 GB peak, 36 min.
- Output: 2,096,155 nuclei × 17,266 genes, 97 donors, 444 donor × region units.

## Source
- File: `$D/multiregion_sc/adata_multiregion_prashant.h5ad`. This is the CELLxGENE release with all nuclei.
- `X` is log-normalised. **Raw integer counts are in `raw/X`**, with the same cells and the same 17,266 genes. `raw/X` is streamed in row chunks to float32.
- The `*_2pct` / `*_10pct` (`_raw`) files are random subsets of the same object and are not used.
- Only 17,266 genes are present, because the CELLxGENE release is gene-filtered. Symbols come from `feature_name`; Ensembl IDs are in `var.gene_ids`. Four duplicate symbols are made unique.

## obs
- `donor` = donor_id; `sample` = donor_region. The hashtag pool is not in the release.
- `region`:

  | tissue | region |
  |---|---|
  | dorsal motor nucleus of vagus | DMNX |
  | medial globus pallidus | GPi |
  | primary motor cortex | motor |
  | prefrontal cortex (BA9) | DLPFC |
  | primary visual cortex | visual |

- `disease`: Parkinson disease → PD; normal → control. Dementia is not given, so no PDD.
- `braak_lb` = `path_braak_lb`.
- `age` comes from development_stage; `pmi` = PMI.
- Extra columns: `brain_bank`, `rin`, `ancestry`, `cell_type_cellxgene`.

## Counts vs paper (claims/notes/prashant2024.md)
| | paper | harmonised |
|---|---|---|
| nuclei after QC | 2,096,155 | 2,096,155 |
| donors with snRNA | 97 (of 100; 75 PD / 25 control overall) | 97: 73 PD / 24 control |
| sex | 62 M / 38 F (100 donors) | 60 M / 37 F |
| donors with all 5 regions | 71 | 71 |

Donors per region × disease:

| region | PD | control |
|---|---|---|
| DLPFC | 70 | 20 |
| DMNX | 70 | 19 |
| GPi | 70 | 19 |
| motor | 70 | 22 |
| visual | 63 | 21 |

- The 3 donors without snRNA are presumably 2 PD and 1 control.
- 5 controls have `braak_lb` 1–2: 4 at stage 1 and 1 at stage 2 (incidental Lewy pathology). They stay `control`.
- PD donors span `braak_lb` 2–6.

Cells per class:

| class | nuclei |
|---|---|
| Oligodendrocytes | 1,093,520 |
| Excitatory neurons | 328,856 |
| Astrocytes | 227,452 |
| Inhibitory neurons | 152,393 |
| Microglia | 118,515 |
| OPCs | 112,585 |
| Neurons (other) | 18,979 |
| Pericytes | 17,865 |
| Ependymal | 12,559 |
| Endothelial | 11,070 |
| T cells | 2,361 |

## Label mapping
`author_cell_type` = `derived_class2`, the author's 12 classes.

| author label | class |
|---|---|
| Oligo | Oligodendrocytes |
| EN | Excitatory neurons |
| IN | Inhibitory neurons |
| Astro | Astrocytes |
| Micro_PVM | Microglia (includes perivascular macrophages) |
| OPC | OPCs |
| Mural | Pericytes |
| Endo | Endothelial |
| Ependymal | Ependymal |
| Adaptive | T cells |
| DMNX_Neu, GPI_Neu | Neurons (other) |

DMNX_Neu and GPI_Neu are the region-specific projection neurons; CELLxGENE calls them 'central nervous system neuron'.

## Label check (original harmonise.label_check)
Cell-weighted agreement is 0.752. Three labels are below 0.5:
- **EN (0.02) and IN (0.06):** the marker majority is "Neurons (other)". This is the same pan-neuronal-score artefact as in the other datasets, and the revised `harmonise.py` and `10_relabel.py` address it.
- **Mural (0.42, majority Fibroblasts):** the class mixes pericytes, SMCs and possibly fibroblasts. It stays as Pericytes, as the authors grouped it.

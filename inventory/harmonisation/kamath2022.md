# kamath2022: harmonisation (GSE178265 human, SNpc + caudate)

Loader: `src/pdrep/build/kamath2022.py` (adapted from `revision/scripts/06a_external_build_and_score.py::load_kamath2022_sn`).
Build: SLURM job 41861172 (300 GB requested, 55 GB peak RSS, 22 min), log `logs/run_py_h_kamath2022_41861172.out`.

## Sources
- Counts: `GSE178265_Homo_matrix.mtx.gz` (5 GB gz; genes × cells, 1.65 × 10⁹ entries) streamed with pyarrow; `GSE178265_Homo_bcd.tsv.gz` (434,340 barcodes, prefix = library = GEO sample title); `GSE178265_Homo_features.tsv.gz` (41,625 hg19 gene symbols; **no Ensembl IDs** in the deposit, so `var["gene_ids"]` is absent).
- Metadata: GEO soft, 97 human samples (tissue, sex, age, PMI, disease, status, cause of death, FACS classification). All barcode prefixes match a GEO title.
- Supplementary tables (`MOESM3_ESM.xlsx`) hold smFISH counts, cluster statistics and DE, with no per-cell labels or Braak stages. The Broad Single Cell Portal metadata (per-cell clusters) is not on disk. → no author labels, `AUTHOR_TO_CLASS = {}`; `braak_lb` NaN.

## Scope (differs from the revision code)
All human libraries are kept: SN NR4A2-sorted ("Pos*") and unsorted DAPI/NR4A2-negative libraries, plus caudate libraries. Extra obs:
- `enrichment`: `sorted` (FACS "Positive" = NR4A2+, or "NeuN") / `unsorted` (FACS "Negative" = DAPI or NR4A2-negative).
- `facs_classification` (raw GEO value), `cause_of_death`, `gsm`.
- Region: SNpc → `SN`; caudate nucleus → `caudate`. Library `pPDCN5730NeuN22119` is annotated "Substantia nigra pars compacta" on GEO but its title is a caudate (CN) library of a caudate-only donor (5730), so region is taken from the title (`caudate`). This is the one NeuN-sorted library (6,315 nuclei).
- Disease from GEO `status`: Ctrl → `control`, PD → `PD`, LBD → `LBD`; `disease_raw` = "disease / status". LBD donor 4775 has cause of death "Parkinson's disease" but status LBD; status is used.
- Donor = 4-digit ID from the title (`pPDsHSrSNxi3873d200429PosA` → 3873; `pPDCN4340DAPIA030419` → 4340).

## Counts vs paper
| | paper | harmonised |
|---|---|---|
| SN controls | 8 donors, 184,673 nuclei | 8 donors, 184,748 nuclei |
| SN PD/LBD | 10 donors, 202,810 nuclei | 10 donors (6 PD + 4 LBD), 202,810 nuclei (135,344 PD + 67,466 LBD) |
| caudate controls | 4 donors, 46,872 nuclei | 4 donors, 46,782 nuclei |
| total | 434,355 | 434,340 |
| DA neurons | 22,048 (15,684 control, 6,364 PD/LBD) | 20,015 marker-labelled (15,270 control, 2,280 PD, 2,465 LBD) |

The deposit is the authors' post-QC nuclei set (total within 15 of the paper). The paper does not give the PD/LBD split in the text; GEO gives 6 PD + 4 LBD. Caudate count 46,782 vs 46,872 is likely a digit transposition in the paper. SN control count is 75 above the paper.
Libraries: SN 51 sorted + 38 unsorted; caudate 1 sorted (NeuN) + 7 unsorted.

**Possible shared donor:** caudate donor 3898 and SN control donor 3298 have identical GEO metadata (male, 79 y, PMI 22 h, prostate cancer). They may be the same person with a typo in one ID. IDs are kept as deposited; treat them as possibly non-independent when caudate and SN controls are analysed together.

Donors: SN PD 1963, 2142, 3873, 3887, 4560, 4568; SN LBD 2544, 2561, 2569, 4775; SN control 3298, 3322, 3345, 3346, 3482, 4956, 5610, 6173; caudate control 3839, 3898, 4340, 5730.

## Cell classes (marker clustering in finalise, 56 Leiden clusters)
| class | SN sorted | SN unsorted | caudate sorted | caudate unsorted |
|---|---|---|---|---|
| Oligodendrocytes | 18,756 | 160,255 | 0 | 10,666 |
| Neurons (other) | 89,935 | 4,009 | 6,307 | 15,599 |
| Astrocytes | 3,555 | 29,691 | 5 | 8,287 |
| Microglia | 7,776 | 25,260 | 0 | 2,100 |
| DA neurons | 19,369 | 646 | 0 | 0 |
| OPCs | 502 | 12,573 | 3 | 2,460 |
| Endothelial | 2,680 | 3,347 | 0 | 1,269 |
| Pericytes | 4,251 | 1,704 | 0 | 10 |
| Fibroblasts | 2,211 | 1,038 | 0 | 76 |

Label check: cell-weighted agreement 0.942. Clusters with agree < 0.6: leiden_12 (3,718, Microglia, 0.53), leiden_14 (782, Microglia, 0.55), leiden_2 (3,103, Neurons (other), 0.56). In all three the per-cell majority class is the same as the cluster label, so they are mixed/low-signal clusters, not mislabels.
Limitations: the NR4A2-sorted libraries are mostly non-DA neurons (89,935 `Neurons (other)`), as in the paper. Caudate neurons are all `Neurons (other)`; there are no `MSNs` (same pan-neuronal scoring issue as xu2023). Sorted and unsorted libraries must not be pooled for composition analyses; use `enrichment`.

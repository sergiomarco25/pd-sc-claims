# adams2024 — harmonisation notes

Loader: `src/pdrep/build/adams2024.py` · SLURM job 41860357 (OK, 280 s, 64 GB) · log `logs/run_py_h_adams2024_41860357.out`

## Sources and counts
- GSE193688: 34 per-library Cell Ranger ARC 1.0.1 `filtered_feature_bc_matrix.h5` (one library per donor). Only the
  `Gene Expression` features are read (ATAC peaks dropped). X = integer UMI counts as provided (ARC cell calling).
- Metadata: GEO soft file (group Young/Aged/PD, diagnosis, age, sex, PMI, RIN, repository) — identical to Supp. Table 1
  (MOESM3). Supp. Table 1 flags 3 donors "Poor RNA Quality, Omitted from downstream analysis" (YC5442, C4634,
  P6320); they are dropped, as in the paper.
- QC (uniform): >= 200 genes and < 10% mitochondrial UMIs per nucleus.
- Region: `midbrain` (paper/GEO source: substantia nigra dissected from midbrain blocks; GEO tissue "Human post-mortem
  midbrain"); `region_raw` keeps "substantia nigra". Extra columns: `age_group` (young/aged/pd), `rin`, `brain_bank`.
- Both young and aged controls have `disease = control`; separate them with `age_group`. braak_lb not available (NaN).

## Cohort vs paper (claims/notes/adams2024.md)
| | paper | here |
|---|---|---|
| donors | 31 = 9 young + 8 aged + 14 PD | 31 = 9 young + 8 aged + 14 PD (midbrain: 17 control, 14 PD) |
| nuclei before QC | 82,735 | 82,735 (sum of the 31 filtered matrices, exact match) |
| nuclei after QC | 69,289 | 75,906 (our QC is milder; the authors also removed doublets/multiplets) |
| neurons | 489 | 2,130 "Neurons (other)" — the authors attribute neuron loss to multiplet removal, which we do not do |

Nuclei per library before/after QC are printed in the log (range 1,223–4,251 after QC).

## Labels
No author cell labels are public: cells labelled by Leiden + marker scores in `finalise()` (author_cell_type = leiden id).
Cells per class: Oligodendrocytes 58,166; OPCs 6,468; Astrocytes 4,460; Microglia 4,361; Neurons (other) 2,130;
T cells 321. No Endothelial cluster (paper: 125 EC — too few to form a cluster). Cell-weighted agreement 0.929.
Clusters with agree < 0.5: leiden_14 (OPCs, 0.48; marker majority also OPCs — mixed OPC/oligo), leiden_19 (T cells,
0.35; small, 321 cells, majority still T cells), leiden_21 (Astrocytes, 0.31, majority Oligodendrocytes; 676 cells,
probably astro–oligo doublets). Kept as labelled.

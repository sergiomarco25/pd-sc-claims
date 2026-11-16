# dehestani2024 — harmonisation notes

Loader: `src/pdrep/build/dehestani2024.py` · conversion `scripts/convert/dehestani2024.R` (SLURM job 41861522, run via `slurm/run_r.sbatch`) · build job 41861783.

## Source and conversion
- GEO GSE272760 provides `GSE272760_all_new_SCT_harmony_res0.1_v1.rds.gz`, a Seurat v4 object.
  - The file is **gzip-compressed twice**. A plain `readRDS` fails with "unknown input format", so the script falls back to `readRDS(gzcon(gzfile()))`.
  - The object holds the assays RNA (v4 `Assay`) and SCT, and the reductions pca, harmony and umap.
- The conversion runs in env `mast_rev`, which has no SeuratObject. `slotNames()` and `@` fail without the class definition, so slots are read with `attr()`.
- Exported: RNA `counts` → `raw/dehestani2024/converted/counts.mtx.gz` (genes × cells), plus `genes.tsv`, `barcodes.tsv` and `meta.tsv`.
  - The counts are integer UMI (0 non-integers): 32,235 genes × 88,876 nuclei.
  - These are the authors' QC-passed nuclei only.
  - The 32,235 genes are the GRCh38-2020-A symbols after the authors' gene filtering.
- Donor metadata comes from Table S1 (`raw/supplements/dehestani2024/13041_2024_1128_MOESM1_ESM.xlsx`): mutation, region, Braak LB, sex, age, brain bank code, PMI (h:mm) and pH.

## Cohort vs paper
| | paper | harmonised |
|---|---|---|
| nuclei | 88,876 after QC | 88,876 |
| donors | 2 LRRK2-PD, 2 GBA1-PD, 2 HC; all male, 65–80 y | 6 donors: LRRK2 ×2, GBA1 ×2 (L444P, E326K), HC ×2; all M; ages 65, 72, 72, 76, 77, 91 |
| regions | PFC + ACC from the same donors | 12 donor × region units (PFC and cingulate, 6 each: 4 PD + 2 control) |
- **Age discrepancy:** one control (NBB 2017-109) is 91 y in Table S1, but the text says the age range is 65–80.
- `donor` is the brain bank code (NBB … or BBN007-26819, the UCL Queen Square LRRK2 case), and `sample` is Sample1–12. The region label ACC is mapped to `cingulate`; the original is kept in `region_author`.
- `disease` is PD for the mutation carriers and control for HC. `genotype` holds GBA1, LRRK2 or none, and `mutation` holds the full variant. `braak_lb` is NaN; Table S1 gives "Braak LB 5-6" for all PD donors, which is kept in `braak_lb_group`.
- PMI: the BBN007 case has "88" in the PMI column, taken as 88 h. That is much longer than the NBB cases (4–6 h). pH is in `brain_ph`.

## Labels
`author_cell_type` = CellType (ODC, ExN1–6, InN1–2, Astro, MG, OPC, Vas); `author_cell_type_broad` = BroadCellType. Vas (505 nuclei) is mapped to Endothelial; the authors call it "vascular".

Cells per class: Oligodendrocytes 44,536 · Excitatory 21,317 · Inhibitory 8,319 · Astrocytes 5,621 · Microglia 4,670 · OPCs 3,908 · Endothelial 505.

Label check (cell-weighted 0.639):
- Glia agree: ODC 0.99, OPC 0.96, Astro 0.85, MG 0.63, Vas 0.61.
- **ExN1–6 (0.00–0.12) and InN1–2 (0.06) are low.** Their marker majority is "Neurons (other)". This is a property of `harmonise.label_check`: the pan-neuronal set (RBFOX3/SNAP25/SYT1/STMN2) outscores the excitatory and inhibitory sets in cortical neurons. The same pattern appears in chia2025, so it is not a labelling error. It was reported to the coordinator. The author labels are kept.

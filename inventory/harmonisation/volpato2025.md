# volpato2025 — harmonisation notes

Loader: `src/pdrep/build/volpato2025.py` · build job 41861085 (SLURM, ~4 min) · outputs in `harmonised/volpato2025*`.

## Source files (raw/volpato2025/, GEO GSE329625)
- `GSE329625_processed_gene_matrix.txt.gz`: genes × nuclei, tab-separated, column names are nucleus IDs with an R `X` prefix (stripped in the loader).
- `GSE329625_full_metadata.csv.gz`: one row per nucleus with `Sample_v2` (donor = library), `Disease`, `Sex`, `CellType`, `CellSubType`, `percent.mt`.

## Counts
- Platform: ICELL8 full-length single-nucleus (SMART-seq type chemistry). Takara Mappa pipeline: STAR to GRCh38.99, then Subread counting.
- The values are **integer read counts**, not UMIs and not TPM. A scan of the first 3,000 genes found 0 non-integers in 68 M values, and finalise's integer check passed on the full matrix. The median count is about 180k–430k reads per nucleus.
- According to GEO data_processing, each gene's count is the sum of two parts: exon-derived counts for all protein-coding genes, plus intron-derived counts for the authors' top 10,000 most variable protein-coding genes. So the matrix has protein-coding genes only (19,803). For the top-10k genes, intronic reads are added on top of the exonic counts.
- Because these are read counts, library sizes are about 100× those of 10x data. Pseudobulk is still summed counts. Downstream tests must not mix this dataset's count scale with UMI datasets without per-dataset normalisation.

## Cohort vs paper (claims/notes/volpato2025.md)
| | paper (v1) | harmonised |
|---|---|---|
| donors | 20 (6 PD B5-6, 3 PD B3-4, 2 iLBD B3-4, 2 iLBD B1-2, 7 control) | 20: 9 PD, 4 iLBD, 7 control, same split by Braak group |
| nuclei | 23,034 after QC | 22,750 (all nuclei present in both GEO matrix and metadata) |
- Region: SN for all donors, from ventral SNpc punches. Each donor is one sample/library, so donor equals sample.
- disease: PD_B5-6 and PD_B3-4 map to PD; ILBD_B1-2 and ILBD_B3-4 map to iLBD; CTR maps to control. `disease_raw` keeps the author string.
- `braak_lb` is NaN because only the group range is known. That range is kept in the extra column `braak_lb_group` (e.g. "5-6").
- age and pmi are NaN. The GEO soft file gives a sex and age for each library, but its pairing of library names with disease conflicts with the per-nucleus metadata. For example, library PD703 is listed as "Control" in the soft file and as PD_B3-4 in the metadata. Sex also disagrees for several libraries. Disease and sex are therefore taken from `full_metadata`, which is consistent with the "PD…" library names, and the soft-file ages are not used. Table S1 (per-donor age) was not available.
- **Open:** donor 13068 (ILBD_B1-2) has only 35 nuclei, so most of its donor × class pseudobulk groups fall below the 10-cell minimum.

## Labels
`author_cell_type` = CellSubType. CellType is kept as `author_cell_type_broad`. The mapping is: DaN_0–3 → DA neurons; ODC_0–2 → Oligodendrocytes; Astrocyte_0–3 → Astrocytes; Microglia_0–2 → Microglia; OPC → OPCs; GABA → Inhibitory neurons; Tcells → T cells.

Cells per class: Oligodendrocytes 11,022 · DA neurons 7,420 · Microglia 1,629 · Astrocytes 1,325 · OPCs 872 · Inhibitory neurons 439 · T cells 43.

Label check (cell-weighted agreement 0.564):
- Glia agree well: ODC 0.78–0.98, OPC 0.93, astrocytes 0.57–0.95, microglia 0.52–0.56, T cells 0.58.
- **DaN_0–3 agree at only 0.01–0.08.** By marker majority, most of these nuclei look like oligodendrocytes or "ependymal". GABA (0.03) looks like "Neurons (other)".
- Investigation:
  - In author DaN, TH is detected in only 8–30% of nuclei. SLC6A3 is detected in 44–89%, but also in 26–42% of glial nuclei. ALDH1A1 and AGTR1 are detected in most nuclei of every class.
  - This points to strong cross-nucleus or ambient signal in this ICELL8 dataset, and to weak exonic DA-marker signal in the DaN clusters.
  - The genes that set DaN apart from all other nuclei (Wilcoxon) are long neuronal genes: PTPRT, ASIC2, CLSTN2, CACNA1E, DLGAP2, MYT1L, GRIN2B, GABBR2. Their signal comes largely from the added intronic counts.
  - So the author DaN clusters are neuronal, but they cannot be confirmed as DA by the short exonic marker list. GABA nuclei have high SNAP25/NEFL/STMN2, so they are clearly neurons.
- Decision: keep the author labels, because the mapping matches the authors' annotation. Treat DA-neuron results from this dataset with caution. Note too that 33% of nuclei are DaN, which is far more than in unsorted SN 10x data.

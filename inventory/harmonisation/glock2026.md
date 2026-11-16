# glock2026 — harmonisation notes

Loader: `src/pdrep/build/glock2026.py` · SLURM jobs 41861081 (first build) and 41861174 (final; label-map fix), 128 GB,
~2 min · log `logs/run_py_h_glock2026_41861174.out`. No claims/notes file exists for this paper. The cohort below comes
from GEO GSE295998 and the authors' cluster table.

## Sources and counts
- 24 per-library `GSM*_<lib>_cellbender_filtered.h5` (CellBender v0.2, default params, on Cell Ranger v6 raw matrices,
  GRCh38 with introns). Read with h5py (`matrix/` group). Values are rounded with `np.rint`. In these files the
  non-integer fraction was 0.0000 in every library, so rounding changes nothing (logged per library).
- `GSE295998_UmapCoordinatesAndClusterAssignments.txt.gz`: 168,270 post-QC nuclei, all scDblFinder singlets.
  Barcodes are `LibraryNN_<bc>`, with NN = position of the library in L1..L8, K1..K8, J1..J8 (checked against the
  Library column). **Cells kept = barcodes present in this table**, so this is the authors' QC and doublet
  removal. Every author barcode was found, and the uniform QC (>= 200 genes, < 10% mito) removed none.
  CellBender called 3,438–37,064 droplets per library; the authors kept 1,951–20,109.
- Metadata from the soft file: subject, sex, age, disease state, cohort (01/02/03 = library batches L/K/J),
  pathology (Lewy body limbic/transitional or diffuse/neocortical; "mild ADNC" in some), and Braak & Braak NFT stage.
  These are stored in `lb_pathology`, `cohort` and `braak_nft`. braak_lb is NaN. PMI is not given.
- donor = Subject. Subject01 has two libraries (L1 in cohort 01 and K1 in cohort 02), so there are 23 donors × 1 region
  and 24 samples. Subject22 (PD) had DBS (noted in the pathology string).
- 58,302 features (full annotation incl. MT-RNR1/2). Duplicated symbols are made unique with suffixes (e.g. RGS5 / RGS5-1).

## Cohort vs paper
GEO says 15 late-onset PD + 8 controls and 168,270 nuclei. Here, amygdala: PD 15, control 8; 168,270 nuclei (exact).

## Labels
author_cell_type = Celltype.l3, falling back to Celltype.l2 where l3 is NA (93 Ciliated_glia cells). Celltype.l1 and
l2 are kept as `author_l1` and `author_l2`. Mapping: Oligo.* → Oligodendrocytes; Ex.* → Excitatory; In.* → Inhibitory;
Astro.* → Astrocytes; OPC.* and COP → OPCs; CG.* and Ciliated_glia → Ependymal; Capillary/Venous/Arterial → Endothelial;
Pericyte and SMC → Pericytes; PV_Fibroblast → Fibroblasts; Tcell → T cells; Bcell → Immune (other); Precursor and
Mix → Neurons (other). Mix has l2 = Excitatory, but its subtype is unclear. MG/PVM.* → Microglia, except MG/PVM.8 and
MG/PVM.10 → Immune (other): these two are perivascular-macrophage-like (CD163+ 60–64%, F13A1+ 23–48%, P2RY12+ 7–12%).
Cells per class: Oligodendrocytes 72,670; Excitatory 29,454; Astrocytes 19,334; Microglia 14,364; OPCs 10,986;
Inhibitory 9,569; Ependymal 4,940; Neurons (other) 1,960; Endothelial 1,752; Immune (other) 1,565; Fibroblasts 640;
T cells 554; Pericytes 482.

Cell-weighted agreement is 0.715. The low value comes almost entirely from **every Ex.* and In.* cluster
(agree 0.00–0.13, marker majority "Neurons (other)")**. This is an artefact of `harmonise.label_check`: the
pan-neuronal marker set always out-scores the excitatory and inhibitory sets in the per-cell argmax, so a correct
ExN/InN label can never "agree". Non-neuronal labels agree well (mostly > 0.75). The other clusters with agree < 0.5
are small: Bcell (0.46, majority Immune (other)), MG/PVM.11 (0.43, majority Microglia), Pericyte (0.39) and SMC (0.42),
whose majority is Pericytes, so these are low-purity rather than wrong labels.

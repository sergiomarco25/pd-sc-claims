# feleke2021 — harmonisation notes

Loader: `src/pdrep/build/feleke2021.py` · final SLURM job 41862511 (OK, 746 s, 64 GB) · log `logs/run_py_h_feleke2021_41862511.out`
Earlier attempts: 41861209 (2,000 MC iterations, cancelled), 41861755 (pending on QOS memory, cancelled), 41861760
(10,000 iterations, retain = knee, cancelled mid-run because samples vanished, see below).

## Sources
- `GSE178146_RAW.tar` was extracted into `raw/feleke2021/extracted/<GSM>_<id>/<id>_raw_feature_bc_matrix/`. It holds
  28 Cell Ranger 3.0.2 **raw (unfiltered)** matrices (GRCh38, Ensembl 93), one library per donor.
- MOESM2 (Supp. Table 1): "Sample information" gives disease group, sex, age at death, disease duration, PMI,
  Braak aSN, Braak tau, Thal, ACG aSN score and LP category. "Single-nucleus RNA-seq" gives the authors' viable nuclei
  per sample, which are used here only for comparison (column `author_n_nuclei`).
- Donor IDs follow the supplement (C036, C048); GEO writes them as C36 and C48.

## Cell calling (no filtered matrices or barcodes are public)
The paper used DropletUtils 1.6.1 emptyDrops with lower = 300 UMI, then kept nuclei with < 5% mito and removed
doublets with DoubletFinder (~7%). DropletUtils is not installed in any env, so emptyDrops is re-implemented in Python
in the loader:
- The ambient profile comes from barcodes with total UMI <= 300, with Good-Turing-like smoothing. The
  Dirichlet-multinomial alpha is fitted by maximum likelihood on those barcodes.
- Each tested barcode gets a DM log-likelihood. P-values come from 10,000 Pólya-urn ambient simulations
  (the DropletUtils default) compared at identical totals. Barcodes are called at BH FDR <= 0.001.
  2,000 iterations is not enough: it floors p at 1/2001, and BH then rejects nothing when the tested set is large.
- **retain = barcodeRanks inflection, not the knee.** The knee and inflection come from a smoothing spline of
  log total on log mid-rank over unique totals > 300, with the top 50 points excluded. The knee was unusable in these
  libraries: it lands at 2,000–19,000 UMI and keeps only 100–3,400 barcodes. In many PD/PDD/DLB libraries, barcodes
  with ~300–2,000 UMI are *more* probable under the ambient model than simulated ambient droplets, because the ambient
  pool is nucleus-like. emptyDrops then calls none of them (sig = 0 in 15 of 28 libraries). With retain = knee these
  samples kept 3–90 nuclei, against 3,000–10,000 in the paper. With retain = inflection, every barcode above the
  inflection is kept and emptyDrops only tests the 300–inflection range.
- Then the uniform QC: >= 200 genes, < 10% mito. Doublets are not removed.

## Nuclei per sample vs paper (Total.Viable.Nuclei)
| donor | here (<10% mito) | here if <5% mito | paper | ratio |
|---|---|---|---|---|
| C036 | 9763 | 9533 | 10493 | 0.93 |
| C048 | 11757 | 11341 | 10653 | 1.10 |
| PD060 | 5759 | 4349 | 4173 | 1.38 |
| PD115 | 5287 | 3162 | 3049 | 1.73 |
| PD163 | 8061 | 5522 | 5228 | 1.54 |
| PD294 | 3895 | 3021 | 3319 | 1.17 |
| PD332 | 5709 | 4435 | 4574 | 1.25 |
| PD341 | 9572 | 9520 | 9770 | 0.98 |
| PD366 | 5291 | 3294 | 3078 | 1.72 |
| PD413 | 7831 | 7744 | 7981 | 0.98 |
| PD415 | 7203 | 6395 | 7277 | 0.99 |
| PD416 | 7569 | 6968 | 6676 | 1.13 |
| PD501 | 6687 | 4616 | 4635 | 1.44 |
| PD523 | 6123 | 4646 | 4489 | 1.36 |
| PD531 | 8912 | 8086 | 7729 | 1.15 |
| PD563 | 9541 | 9337 | 9848 | 0.97 |
| PD566 | 5334 | 3876 | 4287 | 1.24 |
| PD666 | 9991 | 9733 | 9890 | 1.01 |
| PD678 | 5081 | 3632 | 3527 | 1.44 |
| PD683 | 8800 | 7807 | 8383 | 1.05 |
| PD706 | 7850 | 7383 | 7877 | 1.00 |
| PD732 | 5503 | 5366 | 7473 | 0.74 |
| PD747 | 6783 | 6774 | 9803 | 0.69 |
| PDC05 | 14844 | 14634 | 13846 | 1.07 |
| PDC22 | 11893 | 10822 | 10730 | 1.11 |
| PDC34 | 8729 | 8624 | 10010 | 0.87 |
| PDC87 | 6508 | 6232 | 7184 | 0.91 |
| PDC91 | 10074 | 9420 | 9966 | 1.01 |
| **total** | **220,350** | 196,272 | **205,948** | 1.07 |

With the paper's 5% mito cut, the per-sample numbers match closely (most within ±10%). The cell calling therefore
reproduces theirs, and the excess in the uniform build is high-mito nuclei, mostly in low-quality PD/PDD/DLB
libraries. Before DE on these libraries, consider a stricter mito filter or a `pct_counts_mt` covariate.
PD747 is under-called (0.69): its spline inflection sat at 9,206 UMI, so only 2,134 barcodes were retained and
emptyDrops rescued another 4,649. PD732 is also low (0.74).

## Cohort vs paper (claims/notes/feleke2021.md)
Cingulate: 7 control, 7 PD, 7 PDD, 7 DLB, matching the paper. disease = the paper's clinical group. braak_lb = Braak
aSN stage. Extra columns: `braak_tau`, `thal`, `disease_duration`, `acg_asn_score`, `lpc_category`.

## Labels
No public author labels. Marker clustering in the build (job 41862511 ran with the harmonise.py that was current at
submission) gave: Neurons (other) 113,925; Oligodendrocytes 50,189; Astrocytes 27,352; OPCs 12,967; Microglia 11,378;
Endothelial 4,539. Cell-weighted agreement 0.85. The paper's ExN 75,826 / InN 26,467 split is lost, from the same
pan-neuronal issue as jin2024. The coordinator is relabelling with the revised harmonise.py (`scripts/10_relabel.py`).
One cluster has agree < 0.5: leiden_6, labelled Astrocytes (0.24, marker majority Neurons (other); 4,184 cells
checked). It is probably neuron–astrocyte doublets or ambient-contaminated nuclei; check it after relabelling.

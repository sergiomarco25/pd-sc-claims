# feleke2021 — Feleke et al. 2021, Acta Neuropathologica

- DOI 10.1007/s00401-021-02343-x; PMID 34309761; PMCID PMC8357687
- Source used: Europe PMC full-text XML (open access). Abstract, all Results subsections, Discussion, Methods, Fig. 1-6 legends read. Supplementary files (DOCX + 12 XLSX) not opened.

## Cohort as reported
- Region: anterior cingulate cortex (Parkinson's UK Tissue Bank).
- 7 controls, 7 PD (no cognitive impairment), 7 PDD, 7 DLB.
- Sex (female): control 1/7, PD 5/7, PDD 2/7, DLB 0/7 (chi-squared p=0.020) — sex imbalance, especially PD vs control.
- Median disease duration: PD 12 y, PDD 11 y, DLB 6 y (p=0.0099).
- 205,948 nuclei (10x; ~1,398 genes/nucleus) + 24 bulk RNA-seq samples. Cell types: 75,826 ExN, 26,467 InN, 46,662 oligodendrocytes, 25,726 astrocytes, 13,788 microglia, 12,497 OPCs, 4,532 vascular.

## Donor metadata available (Supplementary Table 1)
- Age at death, sex, PMI (mostly <24 h), disease duration, Braak alpha-synuclein stage, Braak tau stage (all <=3 by design), ACC alpha-synuclein semi-quantitative score (0-3, Alafuzoff), Lewy pathology consensus category, RIN.
- Clinical diagnoses from retrospective case notes; 1-year rule for PDD vs DLB.

## Ambiguities / caveats
- snRNA DE: MAST LRT per cell type with age, PMI, sex covariates and cellular detection rate; no donor random effect, so cells are the unit (pseudoreplication risk). |FC|>1.5, FDR<0.05.
- Composition claims use Scaden deconvolution of bulk RNA-seq (trained on the matched snRNA), Wilcoxon on donors — not snRNA nucleus fractions. snRNA vs Scaden proportions correlate poorly overall (rho=0.25) but well for microglia (rho=0.79).
- Microglial increase replicated in an external PD bulk dataset (Scaden).
- Excluded as bulk-only (schema): differential splicing (4,656 intron clusters in 3,751 genes; most in DLB comparisons; DS genes enriched in oligodendrocytes via EWCE using snRNA expression), RBP motif enrichment (GGGGGGG / RBM25 in PDD; HNRNPM motif in DLB; SRSF9 and PCBP2 motifs DLB vs PDD), bulk PC1 pathways. These are headline abstract claims (ii, iii) but are bulk-tissue analyses.
- Many claims compare PDD or DLB rather than PD; recorded in `comparison`.
- SNCA distribution claim (C12) is descriptive without a test.

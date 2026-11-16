# volpato2025: Volpato et al., medRxiv preprint, doi:10.1101/2025.05.28.25328401

**Source used**
- Text: the medRxiv **v1** PDF (posted 2025-05-28), from https://www.medrxiv.org/content/medrxiv/early/2025/05/28/2025.05.28.25328401.full.pdf, extracted with pypdf. The extracted text has broken spacing but is readable.
- Read the abstract, all Results sections, the Methods (composition, DE, CellChat, histology) and the Discussion.
- The medRxiv HTML and JATS XML returned 403, and the Europe PMC PPR1027401 full text was unavailable.
- A **v2** (2025-11-14) exists with a new title ("...overlap with insulin resistance"). Its abstract (from the bioRxiv API) adds "in pODCs, this association reflects insulin resistance with downregulated PI3K-AKT signalling". That analysis is not in v1. It is recorded as C27 with a caveat. Re-check against v2 when it is accessible.

**Cohort as reported (v1)**
- Ventral SNpc, full-length SMART-seq of single nuclei, 2 punches per donor.
- 20 donors in total:
  - 6 idiopathic PD at Braak 5-6
  - 3 PD at Braak 3-4
  - 2 iLBD at Braak 3-4
  - 2 iLBD at Braak 1-2
  - 7 healthy controls
- Mean age 83.
- 23,885 nuclei sequenced and 23,034 passed QC (about 1,194 per donor; about 8,653 genes and 351,886 reads per nucleus).
- Histology (RNAscope/IHC for AGTR1; perineuronal oligodendrocyte counts): 4 PD (Braak 5-6), 4 iLBD (Braak 1-4) and 5 controls.
- The comorbid PD+T2D GWAS used OPDC (965 PD, 69 with diabetes), Tracking Parkinson's (1,436 PD, 97 with diabetes) and UK Biobank (5,460 PD, 866 with diabetes). This is not single-cell work and is not recorded as claims.

**Donor metadata:** Braak stage is given per group. Table S1 has per-donor age and other details; sex, PMI and disease duration are not stated in the main text.

**Ambiguities and caveats**
- Composition: the main test is scProportionTest (a cell-label permutation, i.e. cell-level), intersected with propeller (donor-level).
- DE: Seurat Wilcoxon on cells, which is pseudoreplication. No DE was possible for DaN_3 or Astrocyte_2 because they are "complete[ly] lost/absent" in PD samples.
- "Upregulation of GLP1R/GCGR/GIPR in DaN_3" reads as cell-type specificity, not PD vs control. These genes are low-abundance in 3' data.
- The DTU and isoform-trajectory claims need full-length data and cannot be tested in 10x datasets.
- "AGTR1+ is a state, not a type" is interpretive. The supporting evidence is a pseudotime trajectory plus an iPSC toxin experiment.

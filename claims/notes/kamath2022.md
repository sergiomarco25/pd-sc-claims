# kamath2022: Kamath et al. 2022, Nat Neurosci, doi:10.1038/s41593-022-01061-1, PMID 35513515

**Source used:** Europe PMC full-text XML for PMC9076534 (https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9076534/fullTextXML). Read the abstract, all Results subsections, the main and Extended Data figure legends, and the Methods (MASC, MAST, MAGMA, GSEA). The PubMed MCP text drops italic gene names, so it was used only for cross-checking.

**Cohort as reported**
- SNpc snRNA-seq (10x 3'), with NR4A2-sorted (FANS) and NR4A2-negative libraries:
  - 8 neurotypical controls: 184,673 nuclei; 15,684 DA nuclei used for the taxonomy.
  - 10 PD/LBD donors: 202,810 nuclei, including 6,364 DA neurons.
  - Total 387,483 nuclei and 22,048 DA neurons.
  - The PD/LBD split is in Supplementary Table 2. The text does not give how many donors are PD and how many are LBD.
- Dorsal striatum (caudate) snRNA-seq: 4 neuropathologically normal donors, 46,872 nuclei. Used only for the GWAS enrichment.
- smFISH validation: an extra 10 PD and 10 control frozen midbrains, with 5,339 DA neurons counted.
- Slide-seq: macaque only, 27 arrays.

**Donor metadata in the paper**
- Age and PMI were matched (Extended Data Fig. 7e,f; Welch t-test P = 0.61 and 0.18).
- Sex differs between the arms (Extended Data Fig. 7g). Sex was included as a covariate in MASC and MAST.
- Collection year is reported.
- Midbrain loss scores (0-3) are given where the brain bank reported them (Supplementary Table 2).
- No Braak stages or disease duration in the main text. Supplementary Tables 1-3 may hold more.

**Ambiguities and caveats**
- MASC is a cell-level binomial GLMM with a **library** random effect, not a donor random effect. MAST DE also uses a library random effect.
- Composition was measured within sorted DA neurons. NR4A2 is downregulated in PD, which could bias sorting; the authors addressed this with smFISH (C12, C13).
- 11 of 68 subclusters changed in Fig. 3a. Only the DA, GPNMB microglia and VIM_LHX2 astrocyte clusters are discussed in the text. The OPC and endothelial clusters labelled in Fig. 3a are not recorded as claims.
- The 26 familial PD genes (ref. 38) and the 13 significant TFs are only partly named in the text.
- The cluster name appears both as "CALB_RBP4" and as a CALB1_* name. The full list of 10 DA cluster names is in Fig. 1d.

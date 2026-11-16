# menon2026 — Menon et al. 2026, bioRxiv preprint

Multi-region single-nucleus multi-omic atlas with WGS for rare noncoding variant discovery in PD (Banner Sun Health cohort).
DOI 10.64898/2026.03.05.709922 · PMID 41867822 · PMCID PMC13001318

## Source used
- Europe PMC full-text XML (`/europepmc/webservices/rest/PMC13001318/fullTextXML`), the PMC preprint copy, converted to text.
- Read: abstract, all Results subsections, Discussion, Fig. 1 and Fig. 2 legends, and the Methods on cohort, pooling, DE, composition and LB-stage analysis.
- Fig. 3–6 legends (QTL, ML, rare variants, familial) were skimmed in the Results text only.
- Supplementary tables and figures were not opened.

## Cohort as reported
- Banner Sun Health Research Institute Brain and Body Donation Program (Arizona Study of Aging and Neurodegenerative Disorders).
- Donors: 101 = 80 PD + 21 neurologically normal controls, matched for age and sex. All have 30x WGS.
- Five regions per donor (505 samples):
  - substantia nigra (2 × 50 µm cryosections)
  - putamen
  - anterior cingulate gyrus (BA24/25/32/33)
  - middle temporal gyrus (BA21)
  - cerebellum
- Modality: 10x Multiome. Pooled designs of about 12 samples per batch, genetically demultiplexed.
- Nuclei counts: 4.2 M with good snATAC, 4.5 M with good snRNA, 3.3 M multiome intersection. That is about 32,956 nuclei per individual across regions, and 92% of samples have more than 3,500 nuclei.
- Cell types: 30 subtypes in 7 broad classes.
- Micro-C: FANS-sorted neurons, oligodendrocytes, microglia and astrocytes from MTG, pooled across donors.
- Genetic cohorts (external): GP2 Release 8 (3,611 PD / 350 controls) and UK Biobank WGS (3,070 PD / 30,700 controls).

## Donor metadata
- Supplementary Table 1 has the clinical and neuropathologic annotation.
- Text mentions age at death, sex and PMI as covariates.
- Unified Staging System for Lewy Body Disorders: Stage I N = 21, IIa N = 6, III N = 37, IV N = 36. That sums to 100, and Stage I (olfactory/amygdala LB) has N = 21, the same as the number of controls. Whether Stage I is the controls or an incidental-LB group is ambiguous.
- Disease duration is not mentioned in the text.

## Statistics
- DE: edgeR glmLRT on pseudobulk, with covariates (age, sex, PMI) included only when they correlate with pseudobulk PCs. DEG burden per cell type is checked against 100,000 label permutations.
  - Donors or samples are the unit, so there is no cell-level pseudoreplication.
  - Unclear whether the main DE pools regions per donor ("from each individual within each cell type") or treats donor × region as samples.
- Composition: speckle propeller on per-sample proportions.
- LB stage: edgeR GLM with group + continuous LB stage, LRT on stage. Proportion vs stage uses a linear regression LRT.
- AGTR1 DA subcluster: Wilcoxon on donor proportions within DA neurons.

## Ambiguities / testability
- Most of the paper is genetic (QTLs, ChromBPNet, SKAT rare-variant association, familial pedigrees). These claims (C19, C20) are marked not_testable for sc datasets. Familial-variant nominations (PAQR5, CAMK4) were not recorded as claims.
- The DEG-burden claims (C02–C05, C08–C10) are about how many DEGs a cell type has, not about specific genes. Re-testing means comparing DEG counts under a permutation null in each dataset.
- Specific inhibitory subtypes that increase (C11) and directions for C13 are only in the figures.
- The cell type for the LRFN2 stage decrease (C18) is not stated; MSN is assumed.
- The paper is a preprint (bioRxiv, March 2026).

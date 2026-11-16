# clarence2025 — Clarence et al. 2025, medRxiv preprint

Multi-regional single-nucleus transcriptomic atlas of PD progression (AMP-PD).
DOI 10.1101/2025.07.30.25332436 · PMID 40766144 · PMCID PMC12324657

## Source used
- Europe PMC full-text XML (`/europepmc/webservices/rest/PMC12324657/fullTextXML`), the PMC preprint copy, converted to text.
- Read in full: abstract, all Results subsections, Discussion, Methods (cohort, clustering, crumblr, dreamlet, TPS, LIANA) and all main and Extended Data figure legends.
- Supplementary tables were not opened.

## Cohort as reported
- AMP-PD Release 4. Brain banks: NIH NeuroBioBank (Mount Sinai, Harvard, Miami), University of Miami Brain Endowment Bank and the Udall Center.
- Donor counts conflict:
  - Abstract and Results say 97 donors and 444 samples.
  - Methods say 100 donors = 75 clinically diagnosed PD + 25 neurotypical controls, and also mention "preclinical" cases.
- Donors are aged 57–102, with a male:female ratio of about 3:2 and mainly European ancestry.
- Regions:
  - DMNX (dorsal motor nucleus of vagus)
  - GPi (globus pallidus interna; the MSNs are found here)
  - primary motor cortex
  - prefrontal cortex
  - primary visual cortex
- No SN and no striatum proper. GPi samples contain MSNs, which suggests putamen contamination or adjacent tissue.
- Modality: snRNA-seq only, genotype-demultiplexed with vireo against WGS. 2,092,520 nuclei after QC.
- Taxonomy: 12 classes, 36 subclasses, 62 subtypes. The myeloid compartment has 118,515 nuclei in 13 subtypes.

## Donor metadata
- Braak LB stage 0–6, grouped as control / early (1–2) / medium (3–4) / late (5–6).
- Also available: Braak NFT stage, Hoehn & Yahr, PD medication, visual hallucinations, memory loss, age, sex, PMI, RIN, ancestry, brain bank, and WGS rare-variant carrier status (LRRK2, PRKN, GBA1, PINK1, etc.).
- Disease duration is not mentioned in the text.

## Statistics
- DE: dreamlet pseudobulk LMM.
  - Covariates: age, sex as a random effect, PMI and QC metrics. Brain region is a random effect in the Braak-LB model, but the formula for the PD vs control model omits it.
  - Unit is samples (donor × region). No cell-level pseudoreplication.
- Composition: crumblr with an LMM (brain region and source as random effects; age, RIN and QC as covariates).
- CCI: LIANA rank_aggregate run per donor, then dream LMM on interaction scores, so the test is at donor level.
- TPS: partial correlation of each donor's pseudobulk with the PD DE profile from the same cohort, so it is circular within this dataset.
- Myeloid Braak pseudotime: a neural network predicts Braak LB from cell counts. Its correlation with true Braak LB is modest (r = 0.35 per cell, 0.41 per donor). The module trends are descriptive.

## Ambiguities
- The donor count is 97 or 100 (75 PD / 25 control) depending on the section.
- Directions for subtype composition by stage (C23–C26) are read from text descriptions of Fig. 7c. Exact FDRs are in Supplementary Table 12.
- The GWAS / TF-network "enrichment" (C10) has no formal enrichment statistic in the text, only overlap counts.
- The PsychAD replication of the PD composition effects was not significant (r = 0.51, P = 0.16). This is noted under C29.
- The TPS stage-ranking claims (C14, C15) are descriptive rankings with no formal test.

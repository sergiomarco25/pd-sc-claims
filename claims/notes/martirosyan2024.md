# martirosyan2024: Martirosyan et al. 2024, Mol Neurodegener, doi:10.1186/s13024-023-00699-0, PMID 38245794

**Source used:** Europe PMC full-text XML for PMC10799528. Read the abstract, all Results subsections, the figure legends and the relevant Methods.

**Cohort as reported**
- SNpc snRNA-seq (10x 3' v3, H&E-guided 2 mm punches): 15 sporadic PD and 14 controls, 83,484 nuclei.
- Donors came from the Oregon Brain Bank.
- Spatial: Molecular Cartography (Resolve), 3 PD and 3 controls from the same cohort, with a targeted panel (Supplementary Table 4).

**Donor metadata in the paper**
- Supplementary Table 1 has clinical diagnosis, age, sex, PMI and RIN.
- It also has Lewy body presence in midbrain, limbic (amygdala) and neocortical regions, plus AD changes (neuritic plaques and tangles). Braak is not stated explicitly.
- Most donors are aged 60+, and both sexes are included. Disease duration is not mentioned.

**Ambiguities and caveats**
- Every composition claim uses a **binomial test on nuclei counts** against the overall PD fraction. This treats cells as independent and ignores donors, which is strong pseudoreplication.
- DE used an LMM on cells with a donor random effect (lme4), at BH-adjusted p < 0.1 (a lenient threshold).
- Subpopulation names (Neurons0-5, Astrocytes0-5, Microglia0-5, Oligos0-5) are dataset-specific, so other datasets need label transfer.
- The headline "TH-enriched" glial subpopulations (Astrocytes2, Microglia1, Oligos2) express neuronal and DA genes at low levels. This is likely ambient RNA or doublets from DA neurons, so their "depletion" may only mirror DA neuron loss.
- Oligos3 (OPALIN-high) is reported as increased in PD, the opposite of smajic2022.
- Microglia2 is described as "pro-inflammatory" from its markers alone. This is not recorded as a claim.

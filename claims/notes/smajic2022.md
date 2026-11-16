# smajic2022: Smajić et al. 2022, Brain, doi:10.1093/brain/awab446, PMID 34919646

**Source used:** Europe PMC full-text XML for PMC9050543. Read the abstract, Methods, the full Results and the figure legends. GEO accession: GSE157783.

**Cohort as reported**
- Ventral midbrain snRNA-seq (10x 3' v3.1): 41,435 nuclei in total, 19,002 IPD and 22,433 control.
- **Discrepancy:** the abstract says "six idiopathic PD patients and five controls", but the Results say five IPD cases and six controls. The nuclei counts (more control nuclei) and GEO fit 5 IPD and 6 controls. Check against GEO metadata.
- Dopaminergic neurons are only 0.18% of nuclei. The CADPS2high cluster has 120 cells.
- Immunofluorescence was done on the right hemi-midbrain of the same 11 individuals.
- LCM plus digital PCR validation cohort: C1 and IPD4 from the snRNA cohort, plus 4 new cases and 4 new controls.

**Donor metadata in the paper**
- Mean age: about 77 for IPD (SEM 3) and about 81 for controls (SEM 4).
- PMI: about 22 h for IPD and about 16 h for controls.
- Sex matched (per the abstract).
- Pathology reports describe severe SN neuron loss. Braak stage and disease duration are not in the main text; they may be in Supplementary Table 1.
- All cases had NeuroChip genotyping, with no monogenic cause found.

**Ambiguities and caveats**
- Composition used a t-test on per-sample proportions, at donor level with n = 5 vs 6. The Supplementary Table 6 statistics for astrocytes and oligodendrocytes are not quoted in the text. The Discussion calls the oligodendrocyte decrease in snRNA-seq "a trend".
- DE used the monocle3 quasi-Poisson GLM on cells, with no donor term. This is pseudoreplication.
- Trajectory "enrichment" claims rest on density plots, with no statistic reported.
- The abstract claims glial "proliferation". No Results text supports it (C25), so its testability is uncertain.
- The LCM-dPCR CADPS2 validation is recorded in the notes of C02, not as a separate claim.
- The microglia morphology claim (C06, P = 2e-16) treats cells as independent units.

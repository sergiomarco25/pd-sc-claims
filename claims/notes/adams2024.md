# adams2024 — Adams et al. 2024, Nature Aging

A single-nuclei paired multiomic analysis of the human midbrain reveals age- and Parkinson's disease-associated glial changes.
DOI 10.1038/s43587-024-00583-6 · PMID 38491288 · PMCID PMC11361719

## Source used
- PMC full text (PubMed MCP `get_full_text_article`, PMC11361719) and PMC OAI JATS XML (pmc.ncbi.nlm.nih.gov/api/oai/v1/mh, metadataPrefix=pmc). The XML keeps the gene names (the plain-text API version drops italic gene symbols) and has all main and Extended Data figure legends.
- Supplementary Tables (1–5) were not opened. Gene lists for modules in Methods appear as blanks even in XML for some lists; get them from the publisher PDF if needed.

## Cohort as reported
- Region: substantia nigra, dissected from the midbrain (~3 × 3 × 5 mm blocks). NIH NeuroBioBank (Maryland, Harvard, Miami).
- Modality: 10x Multiome (paired snRNA + snATAC in the same nuclei).
- Donors after QC: 31 = 9 young controls (mean age 24), 8 aged controls (mean 75), 14 PD (mean 81). Three further donors (one per group) were dropped for poor RNA.
- Nuclei: 82,735 before QC, 69,289 after QC. Per type: ODC 45,580; MG 4,482; OPC 4,396; AS 2,428; neurons 489; EC 125; T 110.
- ODC nuclei by group: 15,192 Y / 11,973 A / 18,415 PD. MG: 1,930 / 460 / 2,092. AS: 999 / 397 / 1,032.
- Neurons are strongly under-recovered. The authors attribute this to digitonin permeabilisation causing neuronal multiplets that QC then removed. There are no DA-neuron claims.

## Donor metadata
Supplementary Table 1 is said to hold age, sex and pathology information. The text does not say whether Braak, PMI or disease duration are included. Sex is imbalanced across groups (the authors mention sex genes appearing as DEGs).

## Statistics / pseudoreplication
- Every DE test is a two-sided unpaired t-test on nuclei (FC > 1.25, P < 0.05). No donor-level model is used.
- cPP group comparisons use one-way ANOVA + Tukey on nuclei.
- The composition test (Fig. 1f) is a two-sided t-test. The unit is probably donor proportions, but this is not stated.
- RNA-FISH counts 3–5 fields from each of 3 donors per group, and fields are the test unit.
- The ODC states (healthy / intermediate / disease-associated) come from antimodes of the cPP score distribution. These states are dataset-specific.

## Ambiguities
- The direction of the ODC/MG composition change (aged vs PD) is not given in the text (C12, C13). Fig. 1f or Supplementary Table 2 is needed.
- Several claims are aging-only (young vs aged; C01, C09, C25). They are testable only in datasets with young donors.
- The trajectory claims (C01–C03, C17, C18) rely on a Monocle3 score from combined RNA + ATAC. To replicate them in RNA-only data, use a donor-level signature of the Supplementary Table 4 genes.
- PDE1A (C24) appears only as a FISH "confirmation". The sequencing evidence is presumably in Supplementary Table 5.
- The authors re-applied their cPP pipeline to Smajić et al. 2022 and found a similar disease-associated ODC population (Extended Data Fig. 6c,d). This is recorded in C04 notes.

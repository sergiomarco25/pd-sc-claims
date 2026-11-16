# barbareyes2025 — Oligodendroglia vulnerability in the human dorsal striatum in PD

- DOI 10.1007/s00401-025-02884-5 · PMID 40323467 · PMCID PMC12053221 · Acta Neuropathol 2025
- **Source used:** full text from PMC via the PubMed tool. Gene names in italics were stripped there, so the Europe PMC full-text XML (`/europepmc/webservices/rest/PMC12053221/fullTextXML`) was used for wording. All quotes were checked as verbatim substrings of that text.

## Cohort as reported
- snRNA-seq (10x 3' v3.1; FANS NeuN+/NeuN-, mixed 70/30 or separate lanes): **27 controls** (Pu 27, CN 24), **36 PD** (Pu 35, CN 34). 63 donors: 44 M / 19 F (ctrl 18/9, PD 26/10). Ages: ctrl 25 to >90, PD 60 to >90. Not sex-matched.
- Banks: HBSFRC Los Angeles, Parkinson's UK Brain Bank, MADRC.
- About 200,000 oligodendroglial nuclei (OLs + OPCs), about 3,000 per donor; 15 subpopulations.
- Xenium: **N = 4 donors** (Pu + CN) per Methods. The PD/control split is not stated. 366-gene panel (266 brain panel + 100 custom); 435,157 oligodendroglia cells.
- IHC (contralateral FFPE, MADRC only): TH n = 8 ctrl / 8 PD; α-syn n = 10 ctrl / 8 PD; MBP n not stated.

## Donor metadata available
Supplementary Table 1: age, sex, brain bank and clinical/neuropathological details. The striatum cohort's metadata v2 is in the revision folder (`Human_Samples_Combined_Metadata_v2.xlsx`).

## Ambiguities / caveats
- Composition: scCODA results are shown as pooled percentage differences, with significance given as "p value of confidence < 0.2" (an FDR-like credible-effect threshold). The changes are significant only in putamen; caudate shows a "tendency".
- Most "PD-associated" signature claims compare subpopulations with each other (e.g. PDAO-1 vs MOL-A), not PD vs control within a state. Testing them needs label transfer of the 15 subpopulations or marker-based scoring (Supplementary Table 2).
- The Xenium HSPR-up / myelination-down module claims use a cell-level Wilcoxon test on about 4 donors, so they are effectively anecdotal at donor level.
- The MBP IHC claims (C15–C17) cannot be tested in other single-cell datasets.
- TH loss and α-syn increase in putamen (Fig. 1b) are cohort validation, not extracted as claims.
- PDAO-2's trajectory origin is described as ambiguous and is not extracted.

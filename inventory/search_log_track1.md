# Search log - Track 1 (subcortical: SN/midbrain/VTA, striatum, other subcortical)

Date: 2026-10-21. Searcher: Claude (subagent). Output: search_track1_subcortical.tsv (27 rows: 25 studies + 1 GEO-only dataset + 1 data-resource row).

## PubMed (NCBI E-utilities esearch, all dates up to 2026-10-21)
| # | Query | Hits |
|---|---|---|
| 1 | (Parkinson*[tiab] OR Lewy[tiab]) AND (single-nucleus OR single-nuclei OR single nucleus OR snRNA-seq OR snRNAseq OR single-cell OR scRNA-seq OR multiome OR snATAC OR spatial transcriptom*)[tiab] AND (human OR post-mortem OR postmortem OR patients OR donors)[tiab] | 461 |
| 2 | Parkinson*[tiab] AND (single-nucleus OR snRNA-seq OR single nuclei)[tiab] AND (substantia nigra OR midbrain)[tiab] | 58 |
| 3 | Parkinson*[tiab] AND (single-nucleus OR snRNA-seq OR single-cell)[tiab] AND (putamen OR striatum OR caudate OR accumbens)[tiab] | 65 |
| 4 | Lewy body[tiab] AND (single-nucleus OR snRNA-seq OR single-cell OR spatial transcriptomics)[tiab] | 40 |
| 5 | Parkinson*[tiab] AND (locus coeruleus OR hypothalamus OR thalamus OR globus pallidus OR olfactory bulb OR cerebellum OR vagus OR medulla OR pons)[tiab] AND (single-nucleus OR single-cell OR snRNA-seq)[tiab] | 44 |
| 6 | Parkinson*[tiab] AND (Visium OR Xenium OR MERFISH OR GeoMx OR CosMx OR Slide-seq OR spatial transcriptomic*)[tiab] | 81 |
| 7 | Parkinson*[tiab] AND (snATAC-seq OR single-nucleus ATAC OR chromatin accessibility)[tiab] AND single[tiab] | 17 |
| 8 | Parkinson*[tiab] AND dopaminergic neurons[tiab] AND (single-nucleus OR snRNA-seq)[tiab] AND human[tiab] | 10 |
| 9 | Sarkar[au] AND Parkinson*[tiab] AND (single-nucleus OR single-cell OR snRNA-seq)[tiab] | 2 (none relevant) |
| 10 | Dai[au] AND Parkinson*[tiab] AND (single-nucleus OR snRNA-seq)[tiab] AND nigra[tiab] | 0 |
| 11 | (Parkinson* OR synucleinopath* OR Lewy OR multiple system atrophy)[tiab] AND (single-nucleus OR snRNA-seq OR single nuclei OR single-nuclei)[tiab] AND brain[tiab] | 97 |
| 12 | Parkinson*[tiab] AND (Slide-seq OR Slide-tags OR spatial[ti]) AND (nigra OR midbrain OR putamen OR striatum)[tiab] AND (human OR postmortem OR post-mortem)[tiab] | 14 |
| 13 | Parkinson*[tiab] AND (single-nucleus OR snRNA-seq OR single-cell)[tiab] AND (astrocyte* OR microglia* OR oligodendro* OR vascular OR endothelial OR T cell*)[ti] AND (nigra OR midbrain OR striat* OR putamen)[tiab] | 37 |
| 14 | Chatila Z[au] AND microglia; Volpato C[au] AND Webber C[au]; Graziano M[au] AND Meletis K[au] (publication-status checks) | 8; 0; 3 |
Union of queries 1-8: 583 unique PMIDs; 1-13: 724. Screened by title/abstract (keyword pre-filter + manual review of all post-2018 titles).

## GEO (db=gds, Homo sapiens, GSE)
(Parkinson* OR Lewy) AND (single nucleus OR single-nucleus OR snRNA-seq OR single cell OR single-cell OR spatial OR multiome OR snATAC) -> 106 series. New finds: GSE329625 (Volpato ventral SN), GSE295998 (amygdala PD, Genentech, no paper found). Kamath GSE178265 found by direct lookup (missed by query wording).

## CELLxGENE (curation API, 401 public collections)
Parkinson/Lewy matches: Menon 2026 (16876983), Prashant/N M 2024 (d5d0df8f), Kamath 2022 (b0f0b447), Fullard PFC atlas (cortex, other track), olfactory epithelium Gateway (living tissue, out of scope).

## Web (WebSearch/WebFetch: bioRxiv/medRxiv, ASAP, HCA, Allen)
- Dai/Sarkar SN snRNA (2025) - no match found.
- "Wang et al. Nat Neurosci/Cell 2024-2025 large SN atlas" - only Wang 2024 Sci Adv (GSE184950, Yue lab) found; no larger Wang atlas located.
- ASAP CRN Cloud / Allen ABC Atlas PMDBS release (May 2025): 5 teams (Hafler, Lee, Jakobsson, Scherzer, Hardy), ~3M cells, 211 donors, 9 regions -> resource row.
- bioRxiv: Graziano 2025 putamen (found), Barba-Reyes 2025 preprint (published version used), Chatila 2023 (still preprint), Clarence 2025 medRxiv (Roussos).
- LC / hypothalamus / olfactory bulb / thalamus / cerebellum PD snRNA: none with PD donors found except as regions within multi-region studies (DMNX, GPi in AMP-PD; hypothalamus/SI/VTA in Chatila; cerebellum/putamen in Menon).

## Excluded after screening (selected; reason)
- Weber 2024 eLife (LC snRNA + Visium) - neurotypical donors only, no PD.
- Kilfeather 2024 Cell Rep - mouse.  Khan 2024 PNAS - snRNA in mouse; human data = IHC.  Aguila 2021 - LCM-seq, not single-cell.
- Bae 2023 Mol Brain - reanalysis of Xu 2023 GSE161045.  Wang 2026 ANC, Liu 2025 Curr Med Chem, Fiorini 2024/2025, Gordillo-Gonzalez 2026, Horan-Portelance 2025, Bhattacharya 2025/26, Jung 2026 - reanalyses of public data.
- Kim 2026 bioRxiv (Levin, long-read SV + snRNA) - reuses AMP-PD snRNA of nm2024 (noted in clarence2025 row).
- Cortex-only (other track): Zhu 2024 STM (PFC), Dehestani 2024 (PFC/ACC), Shwab 2024 (temporal cortex), Kang 2025 (cingulate), Chia 2025 (PFC), Liu 2024 eNeuro (PFC), Feleke 2021 (ACC), Goralski 2024, Jin 2026 Acta Neuropathol (temporal cortex), Lee 2026 Nature PsychAD (DLPFC), Langston 2022 (frontal cortex LRRK2 eQTL).
- Hobson 2026 Nat Immunol (leptomeninges + brain immune cells; mostly AD) - not subcortical; flag to other track.
- Stephenson 2025 / Lu 2026 (single-neuron mtDNA, pons) - not transcriptomic.
- Blood/CSF/iPSC/organoid/animal studies - out of scope.

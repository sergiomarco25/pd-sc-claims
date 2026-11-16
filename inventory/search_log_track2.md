# Search log — Track 2: cortical / multi-region / cross-disease PD brain single-cell studies

Date: 2026-10-21. Searcher: Claude (subagent). Output: `search_track2_cortex_multiregion.tsv` (33 rows).

## PubMed (NCBI E-utilities esearch; hit count = total PubMed count; "new" = not seen in earlier queries)
| query | hits | new |
|---|---|---|
| `(Parkinson*[tiab] OR Lewy[tiab] OR synucleinopath*[tiab]) AND (single-nucleus[tiab] OR single-nuclei[tiab] OR snRNA-seq[tiab] OR "single nucleus"[tiab] OR "single-cell RNA"[tiab] OR scRNA-seq[tiab] OR snATAC[tiab] OR multiome[tiab] OR "single-cell transcriptom*"[tiab]) AND (brain[tiab] OR cortex[tiab] OR cortical[tiab] OR postmortem[tiab] OR post-mortem[tiab])` | 254 | 254 |
| `(Parkinson*[tiab] OR Lewy[tiab]) AND (hippocampus[tiab] OR hippocampal[tiab] OR amygdala[tiab] OR entorhinal[tiab] OR "temporal cortex"[tiab] OR "temporal gyrus"[tiab] OR insula[tiab] OR "parietal cortex"[tiab] OR "motor cortex"[tiab] OR "frontal cortex"[tiab] OR cingulate[tiab]) AND (single-nucleus[tiab] OR single-cell[tiab] OR snRNA*[tiab] OR "spatial transcriptom*"[tiab] OR snATAC*[tiab])` | 54 | 33 |
| `(Parkinson*[tiab] OR Lewy[tiab] OR synucleinopath*[tiab]) AND (snATAC[tiab] OR "chromatin accessibility"[tiab] OR multiome[tiab] OR multiomic[tiab]) AND (single-nucleus[tiab] OR single-cell[tiab] OR single-nuclei[tiab]) AND (human[tiab] OR postmortem[tiab] OR post-mortem[tiab])` | 11 | 5 |
| `(Parkinson*[tiab] OR Lewy[tiab]) AND ("spatial transcriptomics"[tiab] OR Visium[tiab] OR Xenium[tiab] OR MERFISH[tiab] OR "GeoMx"[tiab] OR CosMx[tiab]) AND (human[tiab] OR postmortem[tiab])` | 32 | 22 |
| `("dementia with Lewy bodies"[tiab] OR "Lewy body dementia"[tiab] OR "Lewy body disease"[tiab] OR "multiple system atrophy"[tiab] OR "incidental Lewy"[tiab]) AND (single-nucleus[tiab] OR single-cell[tiab] OR snRNA*[tiab] OR single-nuclei[tiab])` | 65 | 26 |
| `Parkinson*[tiab] AND (single-nucleus[tiab] OR single-cell[tiab] OR single-nuclei[tiab]) AND (GWAS[tiab] OR heritability[tiab] OR eQTL[tiab]) AND (brain[tiab] OR cortex[tiab])` | 37 | 9 |
| `Goralski T[au] AND spatial` | 2 | 2 |
| `Chatterjee[au] AND Lewy AND (single OR spatial)` | 2 | 2 |
| `Tuddenham JF[au] AND microglia` | 9 | 9 |
| `Gate D[au] AND Lewy` | 1 | 0 |
| `(Parkinson*[tiab]) AND (single-nucleus[tiab] OR snRNA-seq[tiab]) AND (cortex[tiab] OR cortical[tiab])` | 35 | 0 |
| `(neurodegenerative[tiab] OR "Alzheimer"[tiab]) AND Parkinson*[tiab] AND (single-nucleus[tiab] OR snRNA-seq[tiab]) AND (atlas[tiab] OR cross-disease[tiab] OR "across diseases"[tiab])` | 10 | 1 |

Also via PubMed MCP tool (same day):
- `Parkinson[Title/Abstract] AND (single-nucleus OR single-cell OR snRNA-seq) AND (cortex OR cortical OR cingulate OR prefrontal) AND post-mortem` — 7 hits
- `(Lewy body dementia OR Parkinson's disease dementia OR dementia with Lewy bodies) AND (single-nucleus OR snRNA-seq OR single-cell RNA)` — 90 hits
- `Parkinson AND (single-nucleus OR single-cell) AND (multi-region OR multiregion OR atlas OR "brain regions") AND human brain` — 53 hits

Titles screened manually after removing obvious mouse/iPSC/organoid/blood/CSF/review hits; abstracts fetched (efetch) for ~45 candidates; data-availability statements read from Europe PMC full text where available.

## Europe PMC preprints (SRC:PPR)
- `(Parkinson OR Lewy) AND (single-nucleus OR snRNA-seq OR single-cell OR multiome OR snATAC OR spatial transcriptomics) AND (postmortem OR human brain OR cortex)` — 22 hits; new relevant: Graf 2025 (spatial DLPFC/ACC), Roy 2025 (reanalysis, excluded).

## bioRxiv API (details by DOI) — abstracts verified for
Fairbrother-Browne 2025, Kang 2025, Menon 2026, Breiter 2025, Shwab 2024 preprint, Roy 2025, Graf 2025, Vidyadhara 2024, Pak 2025.

## Web search (WebSearch, extended/standard)
- bioRxiv PD snRNA cortex multi-region atlas 2025
- ASAP CRN PMDBS Team Hardy/Scherzer/Lee data release
- LBD/PDD snRNA cortex bioRxiv/medRxiv 2024–2026
- Menon 2026 PD multi-region atlas CELLxGENE
- large-scale PD cortex snRNA atlas 2025 (Nature/Cell)
- Scherzer PD cortex snRNA temporal cortex; PD5D MTG
- ASAP Team Hafler; Team Hardy (ACC/IPL; 8 regions Braak 3/4); Team Lee
- Toomey/Fairbrother-Browne/Ryten 8 regions; Marshall/Paquette hippocampus-MFG-SN; Columbia Menon/Sims cingulate+SN
- "Molecular and cellular signatures differentiate PD from PDD"
- LBD/PDD snATAC/multiome cortex; Corces MTG PD multiome (IGVF portal)
- GBA/LRRK2 PD cortex snRNA; DLB snRNA hippocampus/entorhinal/amygdala/cingulate
- Tuddenham 2024 living microglia PD donors

## Repositories
- Zenodo API (`q=PMDBS`, `single-nucleus Parkinson postmortem ASAP`, `snRNA Parkinson cortex ASAP`): found ASAP CRN Cloud dataset deposits — 19932738 (Hardy ACG+IPL, 128 samples), 16777109 (Hardy 8 regions, 18 PD Braak 3/4 + 20 ctrl), 19938020/19937390/19937284 (Scherzer PD5D MTG snRNA/Visium), 19937391 (Sulzer/Columbia cingulate+SN), 19936443 (hippocampus/MFG/SN), 20078306 (Aerts multiome ACC/midbrain/M1), 19932739 (Jakobsson SN/PUT/AMY/PFC), 19932730 (Hafler PFC), 18989318/20078073 (PD5D midbrain).
- Allen Brain Cell Atlas ASAP-PMDBS (v20250331): ~3M cells, 9 regions, 211 donors, 5 teams (Hafler, Lee, Jakobsson, Scherzer, Hardy); CC BY 4.0 on AWS (processed matrices) — harmonized collection, not a paper; per-team breakdown not given.
- CELLxGENE: d5d0df8f-… (Prashant 2024 multi-region, 2,096,155 cells); 84ce6837-… (PsychAD DLPFC).
- AMP-PD Release 4 (Prashant/Clarence cohort). Synapse (PsychAD syn60084804; Shwab syn60245188/syn50996869). GEO esummary for GSE202210, GSE235914, GSE303823, GSE288856 (GEO hit rate limit 429 on others).
- IGVF portal: Corces-lab MTG/posterior cingulate/putamen/SN/cerebellum multiome from PD donors (preview, unpublished).

## Excluded after screening (reason)
- Le Bars 2024 Mol Neurobiol — reanalysis of GSE202210 (Zhu) + SEA-AD; no new data.
- Fiorini 2024/2025, Bhattacharya 2025/2026, Chu 2026, Wu 2024 PLoS Biol, Dehkordi 2025, Wang J 2025 (NRR; AMP-PD reanalysis), Gordillo-González 2026, Nakatsuka 2025, Martins-Ferreira 2025 (HuMicA), Palma 2026, Guvenek 2024, Qian 2023, Fan 2023, Alvarado 2023, Pak 2025, Roy 2025/2026, Horan-Portelance 2025 (integrative), Liu Q 2026 (ROSMAP DLPFC), Jung 2026 — reanalyses/meta-analyses of public data.
- Mirzac 2025 (Research; STTT) — FRESH surgical DLPFC from DBS patients (not post-mortem). Flag for caller: out of scope by tissue source but a PD cortex scRNA dataset.
- Gate 2021 Science — scRNA of CSF (brain only histology). Agarwal 2020 — MFG/SN snRNA from controls only. Corces 2020 Nat Genet — scATAC from control brains only (PD in bulk ATAC). Lang 2019 — iPSC. Vidyadhara 2024 — mouse snRNA. Reed 2026 — iPSC. Kon 2023 — no sequencing. Bathe 2024 — GeoMx region-level immune profiling (not single-cell). Kalef-Ezra/Perez-Rodriguez — single-cell DNA CNV. Stephenson 2025 — single-neuron mtDNA. Breiter 2025 — imaging.
- Other tracks (midbrain/striatum): Rydbirk 2026 (striatum MSA/PD), Rumpf 2026 (midbrain spatial), Martirosyan 2024 (SNpc), Smajić 2022, Xu 2023 Nat Aging (striatum), Barba-Reyes 2025/2026 (striatum), Ma 2025 Nat Commun (SN spatial), Shin 2025 (SN spatial), Kamath, Wang Q.
- Seeds not found / unresolved: "Pfaff/Tesla", "Phan", "Gerrits" (Gerrits 2022 is FTD, no PD), "Kim/Garcia multiregion", "Zhang 2025 large-scale PD cortex atlas", "Chatterjee 2023/2024 DLB/PD cortex" — no matching primary PD-brain single-cell paper located. Possible identities: "Menon multiregion" = Menon S 2026 bioRxiv (5 regions, multiome) or Columbia (V. Menon) Team Sulzer cingulate+SN dataset.

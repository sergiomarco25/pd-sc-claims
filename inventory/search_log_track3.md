# Search log — Track 3: non-RNA-only / spatial modalities (human post-mortem PD brain)

Date run: 2026-10-21. Agent: Claude (subagent, track 3). Output: `search_track3_spatial_atac.tsv` (27 rows: 23 included candidates, 1 dataset-only entry, 3 flagged exclude/borderline).

## PubMed (MCP search_articles) queries
| # | Query (abridged) | Hits |
|---|---|---|
| 1 | Parkinson's disease AND (single-nucleus ATAC OR snATAC-seq OR chromatin accessibility single-cell) AND (brain OR substantia nigra OR post-mortem) | 18 |
| 2 | Parkinson's disease AND spatial transcriptomics AND human | 216 (first 60 screened by title; mostly reviews/mouse) |
| 3 | Parkinson AND (multiome OR "single-nucleus multiomic" OR snMultiome) | 6 |
| 4 | (Parkinson OR "Lewy body") AND (Visium OR GeoMx OR Xenium OR MERFISH OR CosMx OR Slide-seq OR Stereo-seq OR "digital spatial profiling" OR "imaging mass cytometry") | 25 |
| 5 | (Parkinson/Lewy body/synucleinopathy)[TiAb] AND (spatial transcriptom*/spatially resolved/spatial profiling/spatial proteomic)[TiAb] AND human/postmortem NOT Review | 50 |
| 6 | (Parkinson OR Lewy body) AND (single-cell methylation / CUT&Tag / long-read single / single-nucleus isoform / CITE-seq / single-cell proteomics) AND brain | 33 |
| 7 | (Parkinson OR Lewy body) AND (chromatin accessibility/snATAC/scATAC/ATAC-seq) AND (nuclei/single-cell) AND human NOT Review | 27 |
| 8 | Title lookup "Characterization of altered molecular mechanisms in PD through cell type-resolved multiomics" | 1 (Lee 2023) |
| 9 | Chatila Z[Author] AND Parkinson | 3 (preprint only; no journal version) |
| 10 | (Parkinson OR Lewy body) AND (SN/midbrain/striatum/putamen/cortex) AND (spatial transcriptomics/Visium/in situ sequencing) AND (snRNA-seq/single-nucleus) AND human, 2021+ | 7 |
| 11 | PD human brain (single-cell/single-nucleus) (spatial/multiome/ATAC/IMC/proteomic) postmortem | 7 |
| 12 | Yoo S[Author] AND Park SH[Author] AND Parkinson AND (single-nucleus OR spatial) | 0 |
| 13 | (PD/LB/synucleinopathy) AND (IMC/CODEX/multiplexed IF/hyperplexed/PhenoCycler/MIBI) AND human brain | 2 |
| 14 | (PD/LB) AND (GeoMx/DSP/LCM single-neuron) AND (transcriptom*/proteom*) AND human | 1 (irrelevant) |
| 15 | (DLB/LBD/MSA/PD) AND (spatial transcriptomic/spatial proteomic/single-nucleus multiome/snATAC/isoform sequencing) AND postmortem NOT Review | 2 |

All titles screened via get_article_metadata; abstracts read for any human-brain candidate.

## Web searches (bioRxiv/medRxiv/Zenodo/ASAP)
- PD snATAC/multiome SN (Lee 2023, Chatila, Adams medRxiv) — found PMC/biorxiv pages.
- Scherzer single-neuron spatial DA neurons — only Kilfeather 2024 (mouse; excluded).
- bioRxiv PD Xenium/MERFISH/CosMx 2025 — CosMx colon PD (excluded, gut); no brain hit.
- Visium SN PD human — Yoo 2024 bioRxiv (new), Ma 2024 bioRxiv (=Ma 2025 Nat Commun), Aguila LCM-seq.
- "Marco Salas Xenium striatum" — Barba-Reyes 2025 Acta Neuropathol; FENS 2026 Cadiz striatal neuron poster (Gonzalez et al., Xenium validation; conference abstract only, not tabulated).
- PD snATAC/multiome 2025–2026 — Menon 2026, Sigalova 2026.
- GeoMx DA neurons PD — Shin 2025 only human hit; others mouse.
- ASAP CRN Cloud spatial — PD5D (Team Scherzer) Visium midbrain Zenodo 10.5281/zenodo.22942068 (100 donors; no paper); Parkinson's Cell Atlas portal (VizIt arXiv 2609.04658) lists Visium + MERFISH datasets (contents not resolvable; PDF not parseable).
- MERFISH/Stereo-seq/Visium HD/CosMx PD human brain — Graf 2025 (Stereo-seq) bioRxiv.
- Martirosyan 2024 — confirmed spatial (Resolve Molecular Cartography, 3 PD/3 ctrl).
- PD cortex Lewy spatial 2025 — Goralski 2024 confirmed; Graf 2025.
- single-nucleus methylation / long-read PD brain — none for methylation; long-read = Liu 2024 eNeuro, Shwab 2025 bioRxiv.
- LBD/PDD snATAC multiome (Hardy/Ryten/ASAP) — none found.
- Xenium PD 2026 preprint — Barba-Reyes 2026 Research Square (own study) via PubMed 42094040.
- 2026 PD spatial preprints — none new beyond the above.

Details (donor counts, accessions) extracted from PMC/Europe PMC full text (OA) where available, else bioRxiv full-text pages.

## Not tabulated (screened out)
Mouse/rat/fly/iPSC/organoid spatial or ATAC (Kilfeather 2024, Gaertner 2025, Horan-Portelance 2024/2026, Jeon 2026, Chen LRRK2-R1441G, Reed 2026 iPSC, Beilina 2026, Hong 2026), gut CosMx PD (2025), Eggers 2025 (LCM proteomics of aging SN, no PD), Dutta 2025 (mouse spatial proteomics), Zhu 2024 STM (snRNA + bulk proteomics; covered by snRNA track), Huuki-Myers 2024 (DLPFC neurotypical), Vicari 2023, imaging-transcriptomics MRI studies, reviews.

## Uncertainties / follow-up
- chatila2023: data availability not located; no journal version found.
- jin2026: LBD cases — whether any are PDD needs Suppl. Table 1.
- bathe2024: donor counts/regions not extracted; LBD group may not be clinical PD.
- chen2021/chen2022/chen2023: control/PD counts partly inferred; chen2022 not OA.
- sigalova2026: per-region PD/control split unverified.
- menon2026: region list inferred from text.
- liu2024: per-group counts not extracted.
- scherzer2026 (PD5D Visium): dataset only, group sizes unknown; watch for paper.
- kamath2022: Slide-seq is control-only; spatial cannot test PD effects.

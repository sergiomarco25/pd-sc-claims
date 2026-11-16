# martirosyan2024 — harmonisation notes

- Loader: `src/pdrep/build/martirosyan2024.py`; build job 41860466 (80 GB, ~3 min).
- Source: `$D/published_datasets/adata_martirosyan.h5ad` (83,484 nuclei × 33,537 genes, Cell Ranger symbols, no Ensembl IDs).
- Counts: `X` is a **dense float32 matrix of raw integer UMI counts** (checked: all non-zero values integer, library sizes 2–12 k). Read in 5,000-row chunks and stored as sparse float32. No layers/raw.
- obs: `donor` = Brain Bank ID, `sample` = Sample ID (one library per donor), `region` = SN for all, `disease` from Clinical diagnosis (Control → control, Parkinson's → PD; disease_raw keeps the string). `braak_lb` = NaN (not reported). Kept as extra columns: `braak_nft`, `cerad`, `lb_midbrain`, `lb_limbic`, `lb_neocortex` (Lewy body presence 0/1, PD donors only; NaN for controls), `rin`, `author_broad` (= obs.cell_type), `brain_bank` = Oregon.

## Counts vs paper
| | paper | harmonised |
|---|---|---|
| nuclei | 83,484 | 83,484 |
| donors | 15 PD + 14 control | 15 PD + 14 control (29 donor × region units) |

Cells per class: Oligodendrocytes 35,032; Astrocytes 20,710; Microglia 12,995; OPCs 6,644; Neurons (other) 3,243; DA neurons 2,615; Endothelial 1,560; T cells 347; Inhibitory neurons 338.

## Label mapping (author_cell_type = obs.SUBTYPES)
- oligosN / astrocytesN / microgliaN, Astro, Micro, OPC, T cells → their classes. The unsubclustered labels (Astro, Micro, Neurons) are cells the authors did not assign to a subpopulation.
- neurons00 (= paper Neurons0, markers TH/SLC6A3/SLC18A2/ALDH1A1 in the supplementary neuron marker table) → DA neurons; neurons03 (GAD1/GAD2/SLC32A1) → Inhibitory neurons; neurons01/02/04/05 and Neurons → Neurons (other). VC (vascular, not split) → Endothelial.
- Note: neurons01 is also TH/SLC6A3/SLC18A2-positive (mean log-norm 0.8/1.1/1.7 vs 2.2/2.6/3.4 in neurons00) with high KCNJ6, so it may be a second DA population; it is kept as Neurons (other) because the paper only calls Neurons0 dopaminergic.

## Label check (harmonised/martirosyan2024_label_check.tsv)
Cell-weighted agreement 0.895. Labels with agree < 0.5:
- neurons00 (DA, 0.12) and neurons03 (Inhibitory, 0.09): marker majority is "Neurons (other)". Manual check of mean log-normalised markers confirms the mapping (neurons00: TH 2.2, SLC6A3 2.6, SLC18A2 3.4; neurons03: GAD1 2.2, GAD2 2.2, SLC32A1 1.3). The disagreement is an artefact of `label_check`: the pan-neuronal "Neurons (other)" score (SNAP25, SYT1, ...) beats the specific neuron-class scores (DDC is near zero in this data). See the harmonise.py note in the build report.
- neurons04 (0.26, majority Oligodendrocytes; MBP/PLP1 high), microglia5 (0.35, Oligodendrocytes) and astrocytes5 (0.48, Oligodendrocytes): small clusters with strong oligodendrocyte signal, likely doublets/ambient. Kept as the authors labelled them.
- VC (0.61) and microglia3 (0.61) are above 0.5.

# smajic2022: harmonisation (GSE157783, ventral midbrain)

Loader: `src/pdrep/build/smajic2022.py` (adapted from `revision/scripts/06a_external_build_and_score.py::load_smajic2022_midbrain`).
Build: SLURM job 41861160 (96 GB, 416 s), log `logs/run_py_h_smajic2022_41861160.out`.

## Sources
- Counts: `GSE157783_IPDCO_hg_midbrain_UMI.tar.gz` (dense genes × cells UMI table, 2.2 GB TSV, read in 2,000-row chunks into sparse). 41,435 nuclei = the authors' post-QC set; 26,737 genes.
- Genes: `_genes.tsv` gives Ensembl IDs in the UMI row order (kept in `var["gene_ids"]`); symbols from the GRCh38 10x features file of GSE161045 (unmapped IDs keep the Ensembl ID).
- Cells: `_cell.tsv` (`barcode`, `cell_ontology`, `patient`), used for `author_cell_type` and donor.
- Donor metadata: GEO soft (diagnosis, sex, age at death, PMI) + Supplementary Table 1 (`raw/supplements/smajic2022/extracted/awab446 Supplementary_data 1.pdf`, p. 7), transcribed into the loader: age at onset, disease duration, Lewy body Braak stage, non-motor features, SN neuron loss (extra obs `age_onset`, `disease_duration`, `nonmotor_features`, `sn_neuron_loss`).

## Donor IDs and the 6/5 vs 5/6 discrepancy
The cell table codes patients `C1..C6` and `PD1..PD5`; GEO titles are `C1..C6`, `IPD1..IPD5`. Nuclei per patient in the cell table match Supplementary Table 1 "Profiled nuclei" exactly (PD1 1,943 = IPD1, PD2 5,925 = IPD2, PD3 3,345, PD4 2,083, PD5 5,706; C1 3,753 … C6 5,758), so `PDk` = `IPDk`; donors use the GEO/table form.
**Resolved:** GEO (11 samples), the cell table and Supplementary Table 1 all give **5 IPD and 6 controls**. The abstract's "six idiopathic PD patients and five controls" is wrong; the Results are right.

## Counts vs paper
| | paper | harmonised |
|---|---|---|
| donors | 5 IPD / 6 control (Results, GEO) | 5 PD-spectrum (4 PD + 1 PDD) / 6 control |
| nuclei | 41,435 (19,002 IPD, 22,433 control) | 41,435 (19,002 / 22,433) |
| DA neurons | 0.18 % (74) | 74 `DaNs` |
| CADPS2high | 120 | 120 `CADPS2+ neurons` |

## Metadata decisions
- Disease: IPD → `PD`, except **IPD3 → `PDD`**, because Supplementary Table 1 lists "Stroke, cognitive impairment, dementia" for IPD3 (HARMONISATION.md rule: PDD when the source says dementia). IPD4 ("cognitive impairment") and IPD5 ("memory impairment") stay `PD`. `disease_raw` keeps the GEO string for all five; analyses of the paper's IPD vs control contrast should pool PD + PDD.
- `braak_lb`: Supplementary Table 1 Lewy body Braak (IPD 5, 6, 6, 6, 5); controls are "-" in the table → NaN (not assumed 0).
- `sample` = donor (one library per donor; barcode suffix `_1.._11` kept in `library_index`).

## Cell classes (author labels)
| author | class | cells | label-check agree |
|---|---|---|---|
| Oligodendrocytes | Oligodendrocytes | 21,268 | 1.00 |
| Astrocytes | Astrocytes | 4,708 | 0.96 |
| Microglia | Microglia | 3,903 | 0.73 |
| Excitatory | Excitatory neurons | 3,037 | 0.01 |
| OPCs | OPCs | 2,754 | 1.00 |
| Endothelial cells | Endothelial | 1,723 | 1.00 |
| Inhibitory | Inhibitory neurons | 1,548 | 0.03 |
| Pericytes | Pericytes | 1,229 | 0.48 |
| GABA | Inhibitory neurons | 535 | 0.00 |
| Ependymal | Ependymal | 536 | 0.91 |
| CADPS2+ neurons | Neurons (other) | 120 | 0.99 |
| DaNs | DA neurons | 74 | 0.32 |

Cell-weighted agreement 0.829. The low agreements are a marker-scoring artefact, not a mapping error: for Excitatory, Inhibitory, GABA and DaNs the per-cell top marker class is `Neurons (other)` (pan-neuronal RBFOX3/SNAP25/SYT1/STMN2 outscores the subtype sets in almost every neuron; see the harmonise.py note in the report). Pericytes (0.48) are split with endothelial/fibroblast markers, as expected for mural cells in snRNA-seq. Author labels are kept.
CADPS2+ neurons are mapped to `Neurons (other)` (the authors describe them as a TH-low, CADPS2-high neuronal population, not DA neurons).

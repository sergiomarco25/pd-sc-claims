# Harmonisation (Phase 3)

Each dataset gets a loader `src/pdrep/build/<dataset>.py` with:
- `AUTHOR_TO_CLASS`: dict author label → one of `pdrep.harmonise.CELL_CLASSES` (`{}` when there are no author labels; cells are then clustered and labelled by markers).
- `load() -> AnnData`: raw integer counts in X (cells × genes, HGNC symbols as var_names; keep Ensembl IDs in `var["gene_ids"]` when available), obs with the columns in `pdrep.harmonise.OBS`.

Build: `sbatch -J h_<dataset> --mem=<GB>G slurm/run_py.sbatch scripts/09_build_dataset.py <dataset>`. Outputs go to
`/lustre/groups/ml01/datasets/projects/2025_pdstriatum_sergio.salas/pd_sc_replication/harmonised/` (see `harmonise.py` docstring).

## obs rules

| column | rule |
|---|---|
| `donor` | the brain donor ID as the authors give it (one donor can have several regions/samples) |
| `sample` | library / sample ID |
| `region` | one of `harmonise.REGIONS` (SN, VTA, midbrain, putamen, caudate, striatum, PFC, DLPFC, cingulate, temporal, MTG, motor, visual, parietal, amygdala, hippocampus, cerebellum, DMNX, GPi, SVZ, other) |
| `disease` | one of control, PD, PDD, DLB, LBD, iLBD, AD, other — use PDD only when the source says dementia; PD donors without that information stay PD |
| `disease_raw` | the source's own diagnosis string |
| `sex` | `M` / `F` / NaN |
| `age`, `pmi` | numbers (years, hours) or NaN |
| `braak_lb` | Lewy body Braak stage 0–6 as a number when available (not tau / NFT Braak) |
| `author_cell_type` | finest author label available |
| extra columns | allowed (e.g. `enrichment` for sorted nuclei, `batch`, `brain_bank`, `braak_nft`, `cell_state_author`) |

Unit of analysis is donor × region (pseudobulk keys donor, region, cell_class).

## Per-dataset notes

Each loader writes nothing else; its author records `inventory/harmonisation/<dataset>.md`: source files, how counts were obtained (filtered matrix, cell calling, conversion from R), donor/cell counts vs the paper (claims/notes/<paper>.md), label mapping decisions, and the label-check agreement.

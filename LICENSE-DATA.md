# Licence for data and text

The code in this repository (`src/`, `scripts/`, `workflow/`, `deploy/`, `tests/*.py`, `config/`) is licensed under the
Apache License 2.0 (`LICENSE`).

Everything else we created is licensed under the **Creative Commons Attribution 4.0 International licence
(CC BY 4.0)**: https://creativecommons.org/licenses/by/4.0/legalcode

That covers:

- the extracted claims, their structure, notes and ratings (`claims/`, including `claims/store/`);
- the claim store and its exports (`claim_store.sqlite`, knowledge graph GraphML/JSON, `tests/fixture/`);
- consensus statements, reference-knowledge records and questions;
- inventory tables (`inventory/`), reports (`reports/`), figures (`figures/`) and the atlas pages.

You may share and adapt them for any purpose, including commercially, if you give appropriate credit (see
`CITATION.cff`), link to the licence and say whether you made changes.

Not covered: quotations from the source papers (short quotes kept for verification, under their publishers' terms),
supplementary tables of other papers, and the single-cell datasets themselves, which keep their original licences and
access conditions (`inventory/datasets.tsv`). Data of datasets listed under `[release] private_datasets` in
`config/project.toml` are not released.

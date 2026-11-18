"""Loader for one dataset: example2024 (replace with a real one; one module per dataset).

Contract (as in the reference project): `load()` returns an AnnData with raw integer counts in X (cells x genes, gene
symbols as var_names) and these obs columns, filling what the source lacks with NaN / "unknown":
  dataset, donor, sample, region, disease, disease_raw, sex, age, pmi, stage, author_cell_type, cell_class
`AUTHOR_TO_CLASS` maps the authors' cell labels to this field's cell classes (config/harmonise.toml cell_classes).
Harmonisation (`python3 -I src/claimstore/cli.py harmonise example2024`) then checks the schema, labels cell classes
and writes the harmonised and pseudobulk files. Low-memory readers for large h5ad files: claimstore.h5io.
"""
AUTHOR_TO_CLASS = {}


def load():
    import anndata as ad  # imported here so the project's checks run without the analysis environment
    import numpy as np
    import pandas as pd
    raise NotImplementedError("read the source files of example2024 and return an AnnData as described above")

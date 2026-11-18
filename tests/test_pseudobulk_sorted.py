"""Pseudobulk export: "ALL" is the tissue, so a dataset of sorted nuclei only has none; mixed datasets keep it."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ad = pytest.importorskip("anndata")
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore import pseudobulk  # noqa: E402


def files(tmp_path, monkeypatch, enrichment, classes):
    n = len(classes)
    donors = [f"d{i // len(set(classes))}" for i in range(n)]
    obs = pd.DataFrame({"donor": donors, "region": "R", "cell_class": classes, "n_cells": 50,
                        "enrichment": enrichment}, index=[str(i) for i in range(n)])
    X = np.random.default_rng(0).poisson(5, (n, 6)).astype(np.float32)
    (tmp_path / "pseudobulk").mkdir()
    ad.AnnData(X, obs=obs).write_h5ad(tmp_path / "pseudobulk" / "ds.h5ad")
    monkeypatch.setattr(pseudobulk, "H", tmp_path)
    monkeypatch.setattr(pseudobulk, "RES", tmp_path / "out")
    monkeypatch.setattr(pseudobulk, "SORTED_CLASSES", ["Micro"])
    monkeypatch.setattr(pseudobulk, "CONTRASTS", ["c"])
    monkeypatch.setattr(pseudobulk, "donor_extras", lambda ds: pd.DataFrame(index=pd.Index(sorted(set(donors)))))
    monkeypatch.setattr(pseudobulk, "unit_groups", lambda meta: pd.DataFrame(
        {"c": ["case", "ref"] * (len(meta) // 2) + ["case"] * (len(meta) % 2)}, index=meta.index))
    pseudobulk.export("ds")
    return sorted(p.name for p in (tmp_path / "out" / "ds" / "input").glob("*.csv.gz"))


def test_sorted_only_has_no_all(tmp_path, monkeypatch):
    assert files(tmp_path, monkeypatch, "sorted", ["Micro"] * 8) == ["R__Micro.csv.gz"]


def test_mixed_keeps_all(tmp_path, monkeypatch):
    classes = ["Micro", "Astro"] * 8
    got = files(tmp_path, monkeypatch, ["sorted", "unsorted"] * 8, classes)
    assert got == ["R__ALL.csv.gz", "R__Astro.csv.gz", "R__Micro.csv.gz"]

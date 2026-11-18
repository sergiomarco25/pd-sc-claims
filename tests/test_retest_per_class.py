"""Signatures applied per cell class ([tests] signatures_per_class): each class gets the rows naming it."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore import retest  # noqa: E402

SIG = pd.DataFrame({"gene": ["A", "B", "C", "D"], "kind": ["de_up", "de_down", "de_up", "named_genes"],
                    "cell_type": ["Neurons", "Neurons", "Glia", ""]})


def calls(monkeypatch, per_class, classes):
    seen = []
    monkeypatch.setattr(retest, "gene_level_test", lambda ds, regs, cls, con, signs: (
        seen.append((tuple(cls), signs)) or pd.DataFrame({"z": [0.0]})))
    retest.gene_level_by_class("d", ["R"], classes, "case_vs_ref", SIG, -1, per_class=per_class)
    return seen


def test_per_class_maps(monkeypatch):
    seen = calls(monkeypatch, True, ["Neurons", "Glia"])
    assert seen == [(("Neurons",), {"A": 1, "B": -1, "D": -1}), (("Glia",), {"C": 1, "D": -1})]


def test_one_map_by_default_or_for_all_cells(monkeypatch):
    assert calls(monkeypatch, False, ["Neurons", "Glia"]) == [(("Neurons", "Glia"), {"A": 1, "B": -1, "C": 1, "D": -1})]
    assert len(calls(monkeypatch, True, ["ALL"])) == 1

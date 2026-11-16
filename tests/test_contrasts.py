"""Contrasts as data (config/contrasts.toml, src/pdrep/contrasts.py; COMMUNITY_PLAN step 2.2).

Needs pandas:  pytest tests/test_contrasts.py
"""
import sys
from pathlib import Path

import pytest

pd = pytest.importorskip("pandas")
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep import config as cfg  # noqa: E402
from pdrep import contrasts  # noqa: E402


def test_pipeline_contrasts_are_defined():
    used = set(cfg.PIPELINE.get("de_contrasts", [])) | set(cfg.PIPELINE.get("gsea_contrasts", []))
    assert used - {"braak"} <= set(contrasts.CONTRASTS)
    assert contrasts.LEWY == {"PD", "PDD", "DLB", "LBD"}


def test_groups_on_example_units():
    m = pd.DataFrame({
        "disease":   ["PD", "PDD", "DLB", "control", "control", "iLBD", "control", "PD"],
        "braak_lb":  [5, 3, None, None, 0, 2, None, None],
        "age_group": ["", "", "", "Aged", "young adult", "", "", ""],
        "genotype":  ["gba N370S", "", "", "", "", "", "", "LRRK2"],
    }, index=[f"u{i}" for i in range(8)])
    g = contrasts.unit_groups(m)
    assert list(g.PDany_vs_control) == ["case", "case", "", "ref", "ref", "", "ref", "case"]
    assert list(g.LBDany_vs_control) == ["case", "case", "case", "ref", "ref", "", "ref", "case"]
    assert list(g.LBlate_vs_control) == ["case", "", "", "ref", "ref", "", "ref", ""]
    assert list(g.LBmid_vs_control) == ["", "case", "", "ref", "ref", "", "ref", ""]
    assert list(g.LBearly_vs_control) == ["", "", "", "ref", "ref", "case", "ref", ""]
    assert list(g.aged_vs_youngcontrol) == ["", "", "", "case", "ref", "", "", ""]
    assert list(g.GBA_vs_control) == ["case", "", "", "ref", "ref", "", "ref", ""]
    assert list(g.LRRK2_vs_control) == ["", "", "", "ref", "ref", "", "ref", "case"]
    assert (g.UPR_vs_other == "").all()  # column absent: nobody in either group
    assert list(g.columns) == list(contrasts.CONTRASTS)


def test_bad_spec_is_rejected(tmp_path):
    bad = [
        '[contrasts.x]\ncase = { disease = "PD" }\n',                                   # no ref
        '[contrasts.x]\ncase = { disease = "PD" }\nref = { disease = { equals = "c" } }\n',  # unknown operator
        '[contrasts.x]\ncase = { lewy_stage = 3 }\nref = { disease = "control" }\n',    # stage needs a range
        '[contrasts.x]\ncase = {}\nref = { disease = "control" }\n',                     # empty side
    ]
    for i, text in enumerate(bad):
        f = tmp_path / f"bad{i}.toml"
        f.write_text(text)
        with pytest.raises(ValueError):
            contrasts.load(f)

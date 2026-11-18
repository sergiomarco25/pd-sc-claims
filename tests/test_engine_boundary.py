"""The engine (src/claimstore) is field-independent: it names no disease, tissue, dataset or paper of a project and
does not import the project's package (claims/COMMUNITY_PLAN.md step 2.1). Facts about a field come from its config.

Standard library only:  python3 -I tests/test_engine_boundary.py   or   pytest tests/
"""
import csv
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from claimstore import config as cfg  # noqa: E402

ENGINE = sorted((ROOT / "src" / "claimstore").glob("*.py"))
ENGINE_R = sorted((ROOT / "src" / "claimstore").glob("*.R"))  # the DE cache
ENGINE_HTML = sorted((ROOT / "src" / "claimstore").glob("*.html"))  # the atlas page
# words that belong to one field (this project's), never to the engine
FIELD = (r"Parkinson|\bPD\b|\bPDD\b|\bDLB\b|Lewy|nigra|dopamin|Braak|striat|putamen|caudate|substantia|synuclein|"
         r"\bSNCA\b|LRRK2|\bGBA\b|Alzheimer|GLIAD|Barba|\bpdrep\b|PDREP|\bbraak|\bMSNs?\b")


def _project_ids():
    ids = set(cfg.DATASETS) | cfg.PRIVATE_DATASETS | cfg.PRIVATE_PAPERS
    papers = ROOT / "inventory" / "papers.tsv"
    if papers.exists():
        with open(papers) as fh:
            ids |= {r["paper_id"] for r in csv.DictReader(fh, delimiter="\t") if r.get("paper_id")}
    return ids


def test_engine_names_nothing_of_a_field():
    pat = re.compile(FIELD + "|" + "|".join(rf"\b{re.escape(i)}\b" for i in sorted(_project_ids())))
    hits = [f"{p.name}:{i}: {m.group(0)}" for p in ENGINE + ENGINE_R + ENGINE_HTML
            for i, line in enumerate(p.read_text().splitlines(), 1) for m in [pat.search(line)] if m]
    assert not hits, f"field-specific names in the engine (move them to the project config): {hits}"


def test_engine_modules_import_only_the_engine():
    for p in ENGINE:
        src = p.read_text()
        assert not re.search(r"^\s*(from|import)\s+pdrep\b", src, re.M), p.name
    assert len(ENGINE) >= 14, [p.name for p in ENGINE]


def test_project_names_still_work():
    import importlib
    for name in ("config", "meta", "store", "kg"):
        assert importlib.import_module(f"pdrep.{name}") is importlib.import_module(f"claimstore.{name}"), name
    # modules that need the analysis environment (anndata, scanpy) are checked by name only
    for old, new in (("pdrep/harmonise.py", "harmonise"), ("pdrep/tests/engine.py", "retest"),
                     ("pdrep/build/_io.py", "h5io")):
        assert f"from claimstore import {new} as _engine" in (ROOT / "src" / old).read_text(), old


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)

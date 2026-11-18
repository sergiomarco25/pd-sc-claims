"""The template makes a working project for another field in one command (scripts/63_new_project.py; step 2.6).

Generates a project into a temporary folder and checks that it holds nothing of this project (no PD terms, dataset or
paper ids, cluster paths), that its claim vocabularies other than regions and cell types are this project's, and that
its own checks pass. Needs pytest, pyyaml, jsonschema, pandas (as CI installs).
"""
import csv
import importlib.util
import re
import sys
import tempfile
import tomllib
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from claimstore import config as cfg  # noqa: E402

pytest.importorskip("jsonschema")
pytest.importorskip("yaml")
pytest.importorskip("pandas")


def _generator():
    spec = importlib.util.spec_from_file_location("new_project", ROOT / "scripts" / "63_new_project.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _this_projects_names():
    ids = set(cfg.DATASETS) | cfg.PRIVATE_DATASETS | cfg.PRIVATE_PAPERS
    with open(ROOT / "inventory" / "papers.tsv") as fh:
        ids |= {r["paper_id"] for r in csv.DictReader(fh, delimiter="\t") if r.get("paper_id")}
    words = (r"Parkinson|\bPD\b|\bPDD\b|Lewy|nigra|dopamin|Braak|striat|putamen|substantia|synuclein|GLIAD|Barba|"
             r"\bpdrep\b|PDREP|/lustre/|/ictstr01/")
    return re.compile(words + "|" + "|".join(rf"\b{re.escape(i)}\b" for i in sorted(ids)))


def test_new_project_is_clean_and_passes_its_checks(tmp_path=None):
    gen = _generator()
    d = Path(tmp_path or tempfile.mkdtemp()) / "ad-sc-claims"
    values = {"name": "ad-sc-claims", "title": "AD single-cell claim store", "description": "test project",
              "package": "ad_sc_claims", "env_prefix": "AD_SC_CLAIMS", "repository": "", "author": "",
              "date": "2026-01-01", "engine_commit": "test", "source": "the reference project"}
    out = gen.generate(d, values)
    for rel in ("config/project.toml", "config/schema.toml", "config/contrasts.toml", "claims/example2024.yaml",
                "claims/schema/claim.schema.json", "claims/RATING.md", "src/claimstore/rating.py",
                "src/claimstore/ENGINE.txt", "src/ad_sc_claims/build/example2024.py", "tests/test_project.py",
                ".github/workflows/ci.yml", "LICENSE", "LICENSE-DATA.md", "CITATION.cff",
                ".claude/skills/claim-replication/SKILL.md",
                ".claude/skills/claim-replication/references/extraction.md"):
        assert (out / rel).is_file(), rel
    assert sorted(p.name for p in (out / "src" / "claimstore").glob("*.py")) == \
        sorted(p.name for p in (ROOT / "src" / "claimstore").glob("*.py"))
    pat = _this_projects_names()
    hits = [f"{p.relative_to(out)}:{i}: {m.group(0)}" for p in out.rglob("*")
            if p.is_file() and "__pycache__" not in p.parts and p.suffix not in (".pyc",)
            for i, line in enumerate(p.read_text().splitlines(), 1) for m in [pat.search(line)] if m]
    assert not hits, hits[:20]
    assert gen.run_checks(out) is True


def test_new_project_builds_rates_and_serves_a_store(tmp_path=None):
    """With the engine alone, a new project merges its claims, maps contexts, builds and rates a store, serves it
    and makes its atlas page."""
    import json
    import os
    import subprocess
    pytest.importorskip("scipy")
    gen = _generator()
    base = Path(tmp_path or tempfile.mkdtemp())
    values = {"name": "new-field", "title": "New field claim store", "description": "test", "package": "new_field",
              "env_prefix": "NEW_FIELD", "repository": "", "author": "", "date": "2026-01-01", "engine_commit": "test",
              "source": "the reference project"}
    out = gen.generate(base / "new-field", values)
    env = {**os.environ, "NEW_FIELD_WORK": str(base / "work")}
    cli = [sys.executable, "-I", str(out / "src" / "claimstore" / "cli.py")]
    for step in ("merge", "contexts", "store", "rating", "store", "atlas"):
        r = subprocess.run(cli + [step], cwd=out, env=env, capture_output=True, text=True)
        assert r.returncode == 0, (step, r.stderr[-800:])
    db = base / "work" / "results" / "store" / "claim_store.sqlite"
    msgs = [{"jsonrpc": "2.0", "id": 1, "method": "initialize", "params": {"protocolVersion": "2025-06-18"}},
            {"jsonrpc": "2.0", "id": 2, "method": "tools/call", "params": {"name": "overview", "arguments": {}}},
            {"jsonrpc": "2.0", "id": 3, "method": "tools/call",
             "params": {"name": "get_claim", "arguments": {"claim_id": "example2024-C01"}}}]
    r = subprocess.run(cli + ["serve"], cwd=out, env={**env, "NEW_FIELD_STORE_DB": str(db)},
                       input="\n".join(json.dumps(m) for m in msgs) + "\n", capture_output=True, text=True, timeout=120)
    replies = {x["id"]: x for x in map(json.loads, r.stdout.splitlines())}
    assert replies[1]["result"]["serverInfo"]["name"] == "new-field"
    for i in (2, 3):
        assert not replies[i]["result"]["isError"], replies[i]["result"]["content"][0]["text"][:300]
    claim = json.loads(replies[3]["result"]["content"][0]["text"])["claim"]
    assert claim["verdict"] == "Untested"
    page = (out / "reports" / "atlas.html").read_text()
    assert "<title>New field claim store</title>" in page and "__TITLE__" not in page and "/*__DATA__*/null" not in page
    data = json.loads(page.split("const D = ", 1)[1].split(";\n", 1)[0])
    assert [c["id"] for c in data["claims"]] == ["example2024-C01"] and data["field"]["title"] == "New field claim store"
    assert (out / "reports" / "atlas.artifact.html").exists()


def test_template_shares_this_projects_claim_vocabularies():
    mine = tomllib.loads((ROOT / "config" / "schema.toml").read_text())
    tpl = tomllib.loads((ROOT / "template" / "config" / "schema.toml").read_text())
    for k in mine["vocab"]:
        if k not in ("region", "cell_type"):
            assert tpl["vocab"][k] == mine["vocab"][k], k
    assert tpl["required"] == mine["required"] and tpl["list_fields"] == mine["list_fields"]
    assert set(tpl["text"]) == set(mine["text"]) and set(tpl["finding_text"]) == set(mine["finding_text"])

"""Release rules and the knowledge graph on the committed test store (COMMUNITY_PLAN steps 1.5 and 3.1).

The fixture is public test data, so it must never contain results from a private dataset (config [release]).
Standard library only:  python3 -I tests/test_release_rules.py   or   pytest tests/
"""
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep import config as cfg  # noqa: E402
from pdrep import kg  # noqa: E402

FIXTURE = ROOT / "tests" / "fixture" / "claim_store.sqlite"


def test_config_has_release_rules():
    assert "barbareyes2026" in cfg.PRIVATE_DATASETS          # GLIAD data private until publication (2026-11-14)
    assert not cfg.is_public("barbareyes2026") and cfg.is_public("xu2023")
    for ds in cfg.PRIVATE_DATASETS:  # in the pipeline here; absent with its loader from the public repository
        assert ds in cfg.DATASETS or not (ROOT / "src" / "pdrep" / "build" / f"{ds}.py").exists()
    assert cfg.PUBLIC_REPOSITORY.startswith("https://github.com/")
    assert "F-axis2-shift" not in cfg.SEPARABLE_FINDINGS  # axis defined with GLIAD donors: re-derived instead
    assert all(f.startswith("F-") for f in cfg.SEPARABLE_FINDINGS)


def test_rederived_findings_are_configured():
    assert "F-axis2-shift" in cfg.REDERIVED_FINDINGS
    assert not set(cfg.REDERIVED_FINDINGS) & cfg.SEPARABLE_FINDINGS
    for fid, spec in cfg.REDERIVED_FINDINGS.items():
        assert fid.startswith("F-") and spec["table"] and not spec["table"].startswith("/")
        assert (ROOT / spec["script"]).is_file(), spec["script"]
    atlas = (ROOT / "scripts" / "51c_store_browser.py").read_text()
    assert cfg.PRIVATE_MODULES and all(f'"{m}"' in atlas for m in cfg.PRIVATE_MODULES)  # names as the programmes use them
    assert {m for m in cfg.PRIVATE_MODULES if "(Barba-Reyes)" in m} == set(cfg.PRIVATE_MODULES)


def test_fixture_holds_no_private_data():
    con = sqlite3.connect(f"file:{FIXTURE}?mode=ro", uri=True)
    for table in ("claim_dataset", "dataset", "donor_group"):
        for ds in cfg.PRIVATE_DATASETS:
            n = con.execute(f"SELECT count(*) FROM {table} WHERE dataset = ?", (ds,)).fetchone()[0]
            assert n == 0, f"{table} has {n} rows of private dataset {ds}"
    for ds in cfg.PRIVATE_DATASETS:
        assert con.execute("SELECT count(*) FROM claim WHERE source_dataset = ?", (ds,)).fetchone()[0] == 0
    assert con.execute("SELECT count(*) FROM claim WHERE claim_id LIKE 'F-%'").fetchone()[0] == 0
    assert con.execute("SELECT count(*) FROM question").fetchone()[0] == 0
    import importlib.util  # and no text names a private dataset, paper or module (the release's own scan)
    spec = importlib.util.spec_from_file_location("release", ROOT / "scripts" / "59_release.py")
    rel = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rel)
    assert rel.private_mentions(FIXTURE) == {}, rel.private_mentions(FIXTURE)


def test_knowledge_graph_builds_without_private_nodes():
    g = kg.build(FIXTURE)
    ids = set(g.nodes)  # node id -> attributes
    assert any(i.startswith("claim:") for i in ids)
    assert not any(i == f"dataset:{ds}" for i in ids for ds in cfg.PRIVATE_DATASETS)



def test_public_release_filter(tmp_path=None):
    """scripts/59_release.py's filter runs on a store and leaves nothing private (ratings are tested separately)."""
    import importlib.util
    import shutil
    import tempfile
    spec = importlib.util.spec_from_file_location("release", ROOT / "scripts" / "59_release.py")
    rel = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rel)
    d = Path(tmp_path or tempfile.mkdtemp())
    db = d / "public.sqlite"
    shutil.copy(FIXTURE, db)
    con = sqlite3.connect(db)
    removed = rel.make_public(con, {"xu2023"}, {"xu2023"})  # pretend xu2023 were private
    assert removed["claims of private papers"] > 0
    assert con.execute("SELECT count(*) FROM claim WHERE paper = 'xu2023'").fetchone()[0] == 0
    assert con.execute("SELECT count(*) FROM claim_dataset WHERE dataset = 'xu2023'").fetchone()[0] == 0
    assert con.execute("SELECT count(*) FROM sqlite_master WHERE name = 'claim_summary'").fetchone()[0] == 0
    rel.apply_ratings(con, {})  # evidence text from public tests only, even without new ratings
    con.close()
    g = kg.build(db)
    assert "dataset:xu2023" not in g.nodes


def test_redact_private_mentions(tmp_path=None):
    """59_release.redact clears private source datasets and drops sentences naming a private dataset."""
    import importlib.util
    import shutil
    import tempfile
    spec = importlib.util.spec_from_file_location("release", ROOT / "scripts" / "59_release.py")
    rel = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(rel)
    d = Path(tmp_path or tempfile.mkdtemp())
    db = d / "public.sqlite"
    shutil.copy(FIXTURE, db)
    con = sqlite3.connect(db)
    cid = con.execute("SELECT claim_id FROM claim LIMIT 1").fetchone()[0]
    con.execute("UPDATE claim SET source_dataset = 'xu2023', notes = ? WHERE claim_id = ?",
                ("Reported for putamen (Fig. 2d). Same cohort as xu2023 (q = 0.04; not after BH). Kept sentence.", cid))
    con.commit()
    n = rel.redact(con, {"xu2023"})  # pretend xu2023 were private
    src, notes = con.execute("SELECT source_dataset, notes FROM claim WHERE claim_id = ?", (cid,)).fetchone()
    assert src is None and notes == "Reported for putamen (Fig. 2d). Kept sentence.", notes
    assert n["claims with a private source dataset cleared"] >= 1 and n["claim text fields redacted"] >= 1
    assert con.execute("SELECT count(*) FROM claim WHERE source_dataset = 'xu2023'").fetchone()[0] == 0
    con.close()


def test_rederive_public_finding(tmp_path=None):
    """scripts/59b_public_findings.py replaces a finding's tests by those of its public table, private datasets out."""
    try:
        import pandas as pd
        import scipy  # noqa: F401
    except ImportError:
        import pytest
        pytest.skip("needs pandas and scipy")
    import importlib.util
    import shutil
    import tempfile
    d = Path(tmp_path or tempfile.mkdtemp())
    db = d / "public.sqlite"
    shutil.copy(FIXTURE, db)
    (d / "sub").mkdir()
    pd.DataFrame({"dataset": ["xu2023", "kamath2022", "smajic2022", "pooled"], "axis": "PC2", "g": [0.9, 0.5, -0.1, None],
                  "se": [0.3, 0.2, 0.4, 0.1], "n_case": [10, 12, 5, None], "n_ctrl": [8, 10, 6, None]}
                 ).to_csv(d / "sub" / "shift.tsv", sep="\t", index=False)
    spec = importlib.util.spec_from_file_location("public_findings", ROOT / "scripts" / "59b_public_findings.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    con = sqlite3.connect(db)
    con.execute("INSERT INTO claim (claim_id, paper, origin, statement) VALUES ('F-test', 'this-study', 'finding', 'x')")
    con.execute("INSERT INTO claim_dataset (claim_id, dataset, role, outcome) VALUES ('F-test', 'old', 'evidence', 'reproduced')")
    f = {"id": "F-test", "direction": "up", "contrast": "PDany_vs_control", "cell_type": "All", "null_claim": False,
         "evidence": {"table": "old.tsv", "filter": {"axis": "PC2", "dataset": ["xu2023", "kamath2022", "smajic2022"]},
                      "columns": {"est": "g", "n_case": "n_case", "n_ref": "n_ctrl"}, "scale": "smd"}}
    saved = m.store.RESULTS
    m.store.RESULTS = d
    try:
        out = m.rederive(con, {"F-test": f}, {"F-test": {"table": "sub/shift.tsv", "script": "s.py"}},
                         is_public=lambda ds: ds != "xu2023")  # pretend xu2023 were private
    finally:
        m.store.RESULTS = saved
    rows = dict(con.execute("SELECT dataset, outcome FROM claim_dataset WHERE claim_id = 'F-test'").fetchall())
    assert rows == {"kamath2022": "reproduced", "smajic2022": "not_reproduced"}, rows
    assert out["F-test"]["tests"] == 2
    assert con.execute("SELECT analysis_results FROM claim WHERE claim_id = 'F-test'").fetchone()[0] == "results/sub/shift.tsv"
    con.close()


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)

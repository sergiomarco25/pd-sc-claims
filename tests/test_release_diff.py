"""Per-claim changelog between releases (scripts/59c_release_diff.py; COMMUNITY_PLAN phase 4 step 5).

Standard library only:  python3 -I tests/test_release_diff.py   or   pytest tests/
"""
import importlib.util
import shutil
import sqlite3
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixture" / "claim_store.sqlite"


def _module():
    spec = importlib.util.spec_from_file_location("release_diff", ROOT / "scripts" / "59c_release_diff.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_diff_finds_each_kind_of_change(tmp_path=None):
    rd = _module()
    d = Path(tmp_path or tempfile.mkdtemp())
    old, new = d / "old.sqlite", d / "new.sqlite"
    shutil.copy(FIXTURE, old)
    shutil.copy(FIXTURE, new)
    con = sqlite3.connect(new)
    ids = [r[0] for r in con.execute("SELECT claim_id FROM claim ORDER BY claim_id LIMIT 3")]
    gone, rerated, edited = ids
    con.execute("DELETE FROM claim WHERE claim_id = ?", (gone,))
    con.execute("UPDATE claim SET verdict = 'Contradicted', certainty = 'Low', est = coalesce(est, 0) + 1 "
                "WHERE claim_id = ?", (rerated,))
    con.execute("UPDATE claim SET statement = statement || ' (edited)' WHERE claim_id = ?", (edited,))
    con.execute("INSERT INTO claim (claim_id, statement, verdict) VALUES ('new2026-C01', 'A new claim.', 'Untested')")
    con.commit()
    con.close()
    rows, labels = rd.diff(old, new)
    got = {(r["claim_id"], r["change"]) for r in rows}
    assert (gone, "removed") in got and ("new2026-C01", "added") in got
    assert (rerated, "verdict") in got and (rerated, "certainty") in got and (rerated, "effect") in got
    assert (edited, "statement") in got
    counts = rd.write(rows, labels, d / "out")
    assert counts["added"] == 1 and counts["removed"] == 1
    text = (d / "out" / "CHANGES.md").read_text()
    assert "## Verdict changed" in text and rerated in text
    assert (d / "out" / "changes.tsv").read_text().splitlines()[0].split("\t") == rd.FIELDS


def test_no_changes_against_itself():
    rows, _ = _module().diff(FIXTURE, FIXTURE)
    assert rows == []


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)

"""The public export (scripts/62_public_repo.py) leaves local tool settings out but keeps the claim-replication skill,
and the skill names nothing private."""
import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
pytest.importorskip("yaml")


def _export():
    spec = importlib.util.spec_from_file_location("public_repo", ROOT / "scripts" / "62_public_repo.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_skill_exported_other_local_settings_left_out():
    m = _export()
    skill = sorted(p.relative_to(ROOT).as_posix() for p in (ROOT / ".claude/skills/claim-replication").rglob("*.md"))
    assert ".claude/skills/claim-replication/SKILL.md" in skill
    assert not any(m.left_out(p) for p in skill)
    for p in (".claude/skills/atlas-questions/SKILL.md", ".claude/settings.local.json", "manuscript/main.md"):
        assert m.left_out(p), p


def test_skill_names_nothing_private():
    m = _export()
    for p in (ROOT / ".claude/skills/claim-replication").rglob("*"):
        if p.is_file():
            text = p.read_text()
            assert not m.PAT.search(text), p
            assert not any(r.search(text) for r in m.MODULES), p
            assert "/lustre/" not in text and "/ictstr01/" not in text, p

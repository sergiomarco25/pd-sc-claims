"""Every claim and finding file validates against its JSON Schema (claims/schema/, COMMUNITY_PLAN step 2.3).

Also checks what JSON Schema cannot: ids match their file, ids are unique, the committed schema files are what
claims/schema/build_schemas.py writes from config/schema.toml, and the merge step uses those same vocabularies.
Needs PyYAML and jsonschema:  pytest tests/test_claim_schema.py   or   python tests/test_claim_schema.py
"""
import importlib.util
import json
from pathlib import Path

import jsonschema
import yaml

ROOT = Path(__file__).resolve().parents[1]
CLAIMS = ROOT / "claims"
SCHEMA_DIR = CLAIMS / "schema"
NOT_PAPERS = {"claims", "findings", "background", "questions"}


def _load(name):
    return json.loads((SCHEMA_DIR / f"{name}.schema.json").read_text())


def _builder():
    spec = importlib.util.spec_from_file_location("build_schemas", SCHEMA_DIR / "build_schemas.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _errors(schema, records, label):
    v = jsonschema.Draft7Validator(schema)
    out = []
    for r in records:
        for e in v.iter_errors(r):
            out.append(f"{label} {r.get('id', '?')}: {'/'.join(map(str, e.path)) or '(record)'}: {e.message[:160]}")
    return out


def test_schema_files_are_current():
    b = _builder()
    assert _load("claim") == b.claim_schema(), "run python3 -I claims/schema/build_schemas.py"
    assert _load("finding") == b.finding_schema(), "run python3 -I claims/schema/build_schemas.py"
    for name in ("claim", "finding"):
        jsonschema.Draft7Validator.check_schema(_load(name))


def test_merge_script_uses_the_schema_vocabularies():
    spec = importlib.util.spec_from_file_location("merge_claims", ROOT / "scripts" / "05_merge_claims.py")
    merge = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(merge)
    b = _builder()
    assert merge.REQUIRED == b.REQUIRED
    for k, allowed in merge.VOCAB.items():
        assert set(allowed) == set(b.VOCAB[k]), k
    assert {"region", "cell_type"} <= set(b.VOCAB) and b.LIST_FIELDS == {"region", "cell_type", "modality"}


def test_paper_claims_validate():
    schema, problems, ids = _load("claim"), [], []
    files = sorted(p for p in CLAIMS.glob("*.yaml") if p.stem not in NOT_PAPERS)
    assert files
    for f in files:
        recs = yaml.safe_load(f.read_text()) or []
        assert isinstance(recs, list), f"{f.name} is not a YAML list"
        problems += _errors(schema, recs, f.name)
        for r in recs:
            if r.get("paper") != f.stem or not str(r.get("id", "")).startswith(f.stem + "-C"):
                problems.append(f"{f.name} {r.get('id')}: id/paper do not match the file name")
            ids.append(r.get("id"))
    dup = sorted({i for i in ids if ids.count(i) > 1})
    assert not dup, f"duplicate ids: {dup}"
    assert not problems, "\n".join(problems[:40])


def test_findings_validate():
    recs = yaml.safe_load((CLAIMS / "findings.yaml").read_text()) or []
    problems = _errors(_load("finding"), recs, "findings.yaml")
    ids = [r.get("id") for r in recs]
    assert len(ids) == len(set(ids)), "duplicate finding ids"
    assert not problems, "\n".join(problems[:40])


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)

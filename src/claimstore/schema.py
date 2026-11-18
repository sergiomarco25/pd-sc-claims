"""Claim schema of a project: JSON Schema (draft-07) for published claims and for the project's own findings, built
from config/schema.toml (required fields, controlled vocabularies, field descriptions), so a new field changes its
vocabularies without code (claims/COMMUNITY_PLAN.md steps 2.3 and 2.6).

  claim.schema.json     one published claim in claims/<paper_id>.yaml
  finding.schema.json   one finding in claims/findings.yaml (same fields, plus origin and evidence)

The schemas' $id points at the project's public repository ([release] public_repository). Claim validation in the
project's merge step uses VOCAB and REQUIRED from here, so the two cannot drift apart.
Standard library only.
"""
import json
import re
import tomllib
from pathlib import Path

from claimstore import config as cfg

SPEC_FILE = cfg.REPO / "config" / "schema.toml"


def load(path=SPEC_FILE):
    with open(path, "rb") as fh:
        return tomllib.load(fh)


SPEC = load() if SPEC_FILE.is_file() else {}
SCHEMA_VERSION = SPEC.get("version", "1.0.0")
REQUIRED = list(SPEC.get("required", []))
LIST_FIELDS = set(SPEC.get("list_fields", []))
VOCAB = {k: list(v) for k, v in SPEC.get("vocab", {}).items()}
TEXT = dict(SPEC.get("text", {}))
FINDING_TEXT = {"id": "F-<short-name>", "paper": "always this-study",
                "test_template": "'finding': estimates come from the finding's own analysis (evidence), not from a "
                                 "re-test template",
                "origin": "record origin", "themes": "themes, ';'-separated", "contrast": "contrast code",
                "null_claim": "true if the finding asserts the absence of an effect",
                "evidence": "where the per-dataset estimates are", "analysis": "how the finding was produced",
                **SPEC.get("finding_text", {})}
BASE = cfg.PUBLIC_REPOSITORY.rstrip("/") + "/blob/main/claims/schema/" if cfg.PUBLIC_REPOSITORY else ""


def list_pattern(values):
    alt = "|".join(re.escape(v) for v in values)
    return rf"^({alt})(\s*[;,]\s*({alt}))*$"


def field(name):
    p = {"description": TEXT[name]}
    if name == "genes":
        p.update(type="array", items={"type": "string", "minLength": 1})
    elif name in LIST_FIELDS:
        p.update(type="string", pattern=list_pattern(VOCAB[name]))
    elif name in VOCAB:
        p.update(type="string", enum=VOCAB[name])
    else:
        p.update(type="string")
    if name in REQUIRED and p.get("type") == "string" and "enum" not in p and "pattern" not in p:
        p["minLength"] = 1
    return p


def claim_schema():
    props = {k: field(k) for k in TEXT}
    props["id"]["pattern"] = r"^[a-z][a-z0-9]*-C[0-9]{2}$"
    return {
        "$schema": "http://json-schema.org/draft-07/schema#",
        "$id": BASE + "claim.schema.json",
        "title": "Published claim",
        "description": f"One conclusion of a published study, as extracted for re-testing (schema {SCHEMA_VERSION}). "
                       "Field meanings: claims/SCHEMA.md.",
        "version": SCHEMA_VERSION,
        "type": "object",
        "required": REQUIRED,
        "properties": props,
        "additionalProperties": False,
    }


def finding_schema():
    s = claim_schema()
    props = s["properties"]
    t = FINDING_TEXT
    props["id"] = {"description": t["id"], "type": "string", "pattern": r"^F-[A-Za-z0-9][A-Za-z0-9-]*$"}
    props["paper"] = {"description": t["paper"], "type": "string", "enum": ["this-study"]}
    props["test_template"] = {"description": t["test_template"], "type": "string",
                              "enum": VOCAB.get("test_template", []) + ["finding"]}
    props.update({
        "origin": {"description": t["origin"], "type": "string", "enum": ["finding"]},
        "themes": {"description": t["themes"], "type": "string"},
        "contrast": {"description": t["contrast"], "type": "string"},
        "null_claim": {"description": t["null_claim"], "type": "boolean"},
        "evidence": {"description": t["evidence"], "type": "object", "required": ["table"],
                     "properties": {"table": {"type": "string", "minLength": 1}}},
        "analysis": {"description": t["analysis"], "type": "object"},
    })
    s.update({"$id": BASE + "finding.schema.json", "title": "Finding of this study",
              "description": f"A result of this study's own analyses, in the claim format (schema {SCHEMA_VERSION}).",
              "required": REQUIRED + ["origin", "evidence"]})
    return s


def write(out_dir=None):
    """Write claim.schema.json and finding.schema.json (default: claims/schema/ of the project)."""
    out_dir = Path(out_dir or cfg.REPO / "claims" / "schema")
    out_dir.mkdir(parents=True, exist_ok=True)
    paths = []
    for name, s in (("claim", claim_schema()), ("finding", finding_schema())):
        out = out_dir / f"{name}.schema.json"
        out.write_text(json.dumps(s, indent=2, ensure_ascii=False) + "\n")
        paths.append(out)
    return paths


def main():
    for p in write():
        print(p)


if __name__ == "__main__":
    main()

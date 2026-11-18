"""Checks of this project (run by CI): the config, the contrasts, the claim schema and every claim file.

Needs pyyaml, jsonschema, pandas:  pytest tests
"""
import json
import sys
from pathlib import Path

import jsonschema
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from claimstore import config as cfg  # noqa: E402
from claimstore import schema  # noqa: E402

NOT_PAPERS = {"claims", "findings", "background", "questions"}


def test_config_describes_the_project():
    about = cfg.describe()
    assert cfg.CONFIG_FILE == ROOT / "config" / "project.toml"
    assert about["name"] and about["server_name"] and cfg.ENV_PREFIX
    assert isinstance(cfg.DATASETS, list) and cfg.PRIVATE_DATASETS <= set(cfg.DATASETS) | cfg.PRIVATE_DATASETS


def test_contrasts_load():
    import pandas as pd
    from claimstore import contrasts
    assert contrasts.CONTRASTS, "define at least one contrast in config/contrasts.toml"
    m = pd.DataFrame({"disease": ["case", "control", "other"]})
    g = contrasts.unit_groups(m)
    assert list(g.columns) == list(contrasts.CONTRASTS)


def harmonise_problems(root):
    """Inconsistencies of config/harmonise.toml with itself and with the stage of config/contrasts.toml."""
    import tomllib
    h = tomllib.loads((root / "config" / "harmonise.toml").read_text())
    st = tomllib.loads((root / "config" / "contrasts.toml").read_text()).get("stage", {})
    lab, classes, markers = h.get("labelling", {}), set(h.get("cell_classes", [])), set(h.get("markers", {}))
    regions, cols = set(h.get("regions", [])), set(h.get("obs", [])) | set(h.get("donor_extras", []))
    subtypes = set(lab.get("neuron_subtypes", []))
    out = [f"obs lacks {c}" for c in ("dataset", "donor", "region", "disease", "author_cell_type", "cell_class")
           if c not in h.get("obs", [])]
    out += [f"cell_classes lacks {c}" for c in sorted(markers | set(h.get("sorted_classes", [])) | {"Unknown"})
            if c not in classes]
    out += [f"neuron subtype {c} has no markers" for c in sorted(subtypes - markers)]
    if subtypes and lab.get("pan_neuronal") not in markers:
        out.append("pan_neuronal must be a class with markers when there are neuron subtypes")
    for key in ("region_only", "require_detect"):
        out += [f"{key}: {c} has no markers" for c in sorted(set(lab.get(key, {})) - markers)]
    out += [f"region_only: {r} not in regions" for v in lab.get("region_only", {}).values() for r in v if r not in regions]
    out += [f"precedence: {a} > {b} not both neuron subtypes" for a, b in lab.get("precedence", {}).items()
            if not {a, b} <= subtypes]
    out += [f"stage column {c} is not in obs or donor_extras" for c in st.get("columns", []) if c not in cols]
    return out


def test_harmonise_config_is_consistent():
    assert not harmonise_problems(ROOT), harmonise_problems(ROOT)


def test_schema_files_are_current():
    for name, built in (("claim", schema.claim_schema()), ("finding", schema.finding_schema())):
        f = ROOT / "claims" / "schema" / f"{name}.schema.json"
        assert json.loads(f.read_text()) == built, "run python3 -I claims/schema/build_schemas.py"
        jsonschema.Draft7Validator.check_schema(built)


def test_claim_files_validate():
    v = jsonschema.Draft7Validator(schema.claim_schema())
    problems, ids = [], []
    for f in sorted(p for p in (ROOT / "claims").glob("*.yaml") if p.stem not in NOT_PAPERS):
        recs = yaml.safe_load(f.read_text()) or []
        for r in recs:
            problems += [f"{f.name} {r.get('id')}: {e.message[:120]}" for e in v.iter_errors(r)]
            if r.get("paper") != f.stem or not str(r.get("id", "")).startswith(f.stem + "-C"):
                problems.append(f"{f.name} {r.get('id')}: id/paper do not match the file name")
            ids.append(r.get("id"))
    assert len(ids) == len(set(ids)), "duplicate claim ids"
    assert not problems, "\n".join(problems[:30])


def test_findings_validate():
    v = jsonschema.Draft7Validator(schema.finding_schema())
    recs = yaml.safe_load((ROOT / "claims" / "findings.yaml").read_text()) or []
    assert not [e.message for r in recs for e in v.iter_errors(r)]

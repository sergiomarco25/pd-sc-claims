"""config/harmonise.toml (claimstore.harmonise: harmonised data schema, cell classes, marker rules) is consistent with
itself and with the stage of config/contrasts.toml, and holds the values the PD pipeline was built with.

Standard library only:  python3 -I tests/test_harmonise_config.py   or   pytest tests/
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


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


def test_harmonise_config_keeps_the_pd_schema():
    import tomllib
    h = tomllib.loads((ROOT / "config" / "harmonise.toml").read_text())
    assert h["obs"] == ["dataset", "donor", "sample", "region", "disease", "disease_raw", "sex", "age", "pmi",
                        "braak_lb", "author_cell_type", "cell_class"]
    assert len(h["cell_classes"]) == 16 and list(h["markers"]) == h["cell_classes"][:15]
    assert h["labelling"]["require_detect"]["DA neurons"] == {"genes": ["TH", "SLC6A3"], "min_fraction": 0.30}


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)

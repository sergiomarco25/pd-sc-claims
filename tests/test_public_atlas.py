"""Public atlas: the audit finds private names in a page's data; the builder and template support the public mode.

Standard library only:  python3 -I tests/test_public_atlas.py   or   pytest tests/
"""
import importlib.util
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep import config as cfg  # noqa: E402


def _audit():
    spec = importlib.util.spec_from_file_location("audit", ROOT / "scripts" / "61_audit_public.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def test_audit_flags_private_names_only():
    a = _audit()
    ds = sorted(cfg.PRIVATE_DATASETS)[0]
    module = sorted(cfg.PRIVATE_MODULES)[0]
    data = {"public": True, "claims": [{"id": f"{ds}-C01", "ds": "xu2023"}, {"id": "barbareyes2025-C01", "ds": "x"}],
            "prog": [{"programme": module, "datasets": "xu2023;kamath2022"}, {"programme": "heat response"}]}
    html = "<script>const D = " + json.dumps(data) + ";\nconst X = 1;</script>"
    got = a.page_data(html)
    assert got == data
    hits = a.find(got, a.patterns())
    paths = sorted({p for _, p, _ in hits})
    assert paths == ["D.claims[].id", "D.prog[].programme"], paths
    assert not any("barbareyes2025" in t for _, _, t in hits)


def test_template_and_builder_have_public_mode():
    tpl = (ROOT / "scripts" / "atlas_template.html").read_text()
    assert "/*__DATA__*/null" in tpl and "D.public" in tpl and "D.links?.question_form" in tpl
    b = (ROOT / "scripts" / "51c_store_browser.py").read_text()
    assert "--public" in b and "pd_atlas_public.html" in b and "PRIVATE_MODULES" in b
    for name in ("22_composition_meta.py", "25_genetics_meta.py", "27_programmes.py"):
        assert "\ndef pool(" in (ROOT / "scripts" / name).read_text(), name
    for name in ("23_de_meta.py", "24_pathways.py"):
        assert '"--public" in sys.argv' in (ROOT / "scripts" / name).read_text(), name


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)

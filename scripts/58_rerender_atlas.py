"""Re-render the atlas with the current template, re-using the data already embedded in reports/pd_atlas.html.

A fast path for template-only changes (layout, 3D view, wording), and for refreshing the knowledge-graph part from the
store, without re-running scripts/51c_store_browser.py (which needs the analysis environment and the biology results).
Writes reports/pd_atlas.html and reports/pd_atlas.artifact.html exactly as 51c does.

Usage: python3 -I scripts/58_rerender_atlas.py [--kg] [--db path/to/claim_store.sqlite] [--public]
  --kg      rebuild the embedded knowledge graph (D.kg) from the store with pdrep.kg, as 51c does (the private atlas)
  --public  re-render the public atlas, reports/pd_atlas_public.html, in place (template changes only)
Standard library only.
"""
import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep import config as cfg  # noqa: E402

MARK = "const D = "
SLOT = "/*__DATA__*/null"


def embedded_data(html: str) -> dict:
    a = html.index(MARK) + len(MARK)
    b = html.index("\n", a)
    return json.loads(html[a:b].rstrip().rstrip(";"))


def render(tpl: str, data: dict) -> tuple[str, str]:
    html = tpl.replace(SLOT, json.dumps(data, separators=(",", ":"), default=str))
    head, body = html.split("<body>", 1)
    art = head[head.index("<title>"):head.index("</head>")] + body.rsplit("</body>", 1)[0]
    return html, art


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--kg", action="store_true", help="rebuild the embedded knowledge graph from the store")
    ap.add_argument("--db", default=str(cfg.STORE_DB), help="claim store (default from config/project.toml)")
    ap.add_argument("--src", help="atlas holding the data (default: reports/pd_atlas.html, or the public atlas)")
    ap.add_argument("--public", action="store_true", help="re-render reports/pd_atlas_public.html in place")
    a = ap.parse_args()
    if a.public and a.kg:
        sys.exit("--kg rebuilds from the full store; rebuild the public atlas with 51c_store_browser.py --public")
    src = Path(a.src) if a.src else ROOT / "reports" / ("pd_atlas_public.html" if a.public else "pd_atlas.html")

    tpl = (ROOT / "scripts" / "atlas_template.html").read_text()
    assert tpl.count(SLOT) == 1, "template must contain the data slot exactly once"
    data = embedded_data(src.read_text())
    if bool(data.get("public")) != a.public:
        sys.exit(f"{src} is {'a public' if data.get('public') else 'the private'} atlas; use --public only for the public one")
    if a.kg:
        from pdrep import kg
        data["kg"] = kg.build(a.db).compact(kg.similarity_links(a.db), kg.topic_clusters(a.db))
    html, art = render(tpl, data)
    out = ROOT / "reports" / ("pd_atlas_public.html" if a.public else "pd_atlas.html")
    out.write_text(html)
    if not a.public:
        (ROOT / "reports" / "pd_atlas.artifact.html").write_text(art)
    print(f"{out} ({out.stat().st_size / 1e6:.1f} MB) re-rendered" + (" with a rebuilt knowledge graph" if a.kg else ""))


if __name__ == "__main__":
    main()

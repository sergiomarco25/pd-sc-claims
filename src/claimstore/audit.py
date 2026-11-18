"""Audit a public atlas page for private material (claims/COMMUNITY_PLAN.md step 5.1).

Usage: python3 -I scripts/61_audit_public.py reports/pd_atlas_public.html
Reads the data embedded in the page and lists every string that names a private dataset or a private paper (config
[release] private_datasets, private_papers) or one of the private papers' gene modules (private_modules), with where
it sits in the data. Exits 1 if there is any. Numbers computed from private data cannot be found this way; they are
kept out by construction (scripts/59_release.py, scripts/60_public_biology.py, 23/24/56 --public).
Standard library only.
"""
import json
import re
import sys
from pathlib import Path

from claimstore import config as cfg

ROOT = cfg.REPO


def page_data(html):
    """The JSON object the template receives in place of /*__DATA__*/null."""
    m = re.search(r"const D = ", html)
    if not m:
        raise ValueError("no embedded data (const D = ...)")
    data, _ = json.JSONDecoder().raw_decode(html, m.end())
    if data is None:
        raise ValueError("the page holds the template's empty data")
    return data


def patterns():
    ids = sorted(cfg.PRIVATE_DATASETS | cfg.PRIVATE_PAPERS)
    pats = [(f"private dataset or paper {i}", re.compile(rf"\b{re.escape(i)}\b")) for i in ids]
    pats += [(f"private module {m}", re.compile(re.escape(m), re.I)) for m in sorted(cfg.PRIVATE_MODULES)]
    return pats


def find(data, pats):
    hits = []

    def walk(x, path):
        if isinstance(x, dict):
            for k, v in x.items():
                walk(k, path + "{key}")
                walk(v, f"{path}.{k}")
        elif isinstance(x, list):
            for v in x:
                walk(v, path + "[]")
        elif isinstance(x, str):
            for name, p in pats:
                if p.search(x):
                    hits.append((name, path, x[:120]))
    walk(data, "D")
    return hits


def main():
    if len(sys.argv) < 2:
        sys.exit("usage: audit PAGE (a public atlas page)")
    page = Path(sys.argv[1])
    data = page_data(page.read_text())
    if not data.get("public"):
        sys.exit(f"{page} is not a public build (no public flag)")
    hits = find(data, patterns())
    by = {}
    for name, path, text in hits:
        by.setdefault((name, path), []).append(text)
    for (name, path), texts in sorted(by.items()):
        print(f"{len(texts):5d}  {name:45s} {path}   e.g. {texts[0]!r}")
    print(f"{page}: {len(hits)} private mentions in {len(by)} places")
    sys.exit(1 if hits else 0)


if __name__ == "__main__":
    main()

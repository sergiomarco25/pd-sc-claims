"""Linked-data exports of a release (claims/COMMUNITY_PLAN.md step 2.5): claims as JSON-LD and the folder as an RO-Crate.

Usage: python3 -I scripts/59d_linked_data.py RELEASE_DIR   (scripts/59_release.py runs it on every release)
  claims.jsonld           every claim as a schema.org ClaimReview: the claim (text, the paper it appears in by DOI,
                          what it is about: regions as UBERON and cell types as Cell Ontology terms, genes by symbol)
                          and our rating as two Ratings (verdict; certainty 1-4 = very low to high), with the evidence
                          summary (cohorts, pooled effect, concerns) as the review text
  ro-crate-metadata.json  RO-Crate 1.1 description of the release folder: licence, authors (CITATION.cff), every file
                          with format and size (checksums: SHA256SUMS, listed in the crate)
Only schema.org terms and ids relative to the release are used, so nothing depends on a namespace we do not own;
nanopublications need persistent IRIs and wait for them. Standard library only.
"""
import argparse
import json
import re
import sqlite3
from pathlib import Path

from claimstore import config as cfg

ROOT = cfg.REPO
CC_BY = "https://creativecommons.org/licenses/by/4.0/"
FORMATS = {".sqlite": "application/vnd.sqlite3", ".tsv": "text/tab-separated-values", ".json": "application/json",
           ".jsonld": "application/ld+json", ".graphml": "application/graphml+xml", ".md": "text/markdown",
           ".cff": "application/x-yaml"}
LATER = ("release.json", "SHA256SUMS")  # written by 59_release.py after this script


def _split(s):
    return [x.strip() for x in str(s or "").split(";") if x.strip()]


def _obo(curie):
    return f"http://purl.obolibrary.org/obo/{curie.replace(':', '_')}" if re.fullmatch(r"[A-Za-z]+:\d+", curie) else None


CERTAINTY = {"Very low": 1, "Low": 2, "Moderate": 3, "High": 4}


def _set(name):
    return {"@type": "DefinedTermSet", "name": name}


def _terms(names, ids, scheme):
    out = []
    names, ids = _split(names), _split(ids)
    for i, name in enumerate(names):
        t = {"@type": "DefinedTerm", "name": name, "inDefinedTermSet": _set(scheme)}
        curie = ids[i] if i < len(ids) else ""
        if _obo(curie):
            t["@id"], t["termCode"] = _obo(curie), curie
        out.append(t)
    return out


def claims_jsonld(db):
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    con.row_factory = sqlite3.Row
    have = {r[1] for r in con.execute("PRAGMA table_info(claim)")}
    src = {r["paper_id"]: r for r in con.execute("SELECT paper_id, title, doi, year, journal_or_preprint FROM source")}
    graph = []
    for c in con.execute("SELECT * FROM claim ORDER BY claim_id"):
        cid = c["claim_id"]
        claim = {"@id": f"#claim-{cid}", "@type": "Claim", "identifier": cid, "text": c["statement"]}
        s = src.get(c["paper"])
        if c["origin"] == "finding":
            claim["appearance"] = {"@id": "./"}  # a finding of this release's own analyses
        elif s is not None and s["doi"]:
            claim["appearance"] = {"@id": f"https://doi.org/{s['doi']}", "@type": "ScholarlyArticle",
                                   "name": s["title"], "datePublished": str(s["year"] or ""),
                                   "isPartOf": ({"@type": "Periodical", "name": s["journal_or_preprint"]}
                                                if s["journal_or_preprint"] else None)}
        about = (_terms(c["regions"], c["region_ids"] if "region_ids" in have else "", "UBERON")
                 + _terms(c["cell_types"], c["cell_type_ids"] if "cell_type_ids" in have else "", "Cell Ontology")
                 + [{"@type": "DefinedTerm", "name": g, "inDefinedTermSet": _set("HGNC symbol")} for g in _split(c["genes"])])
        if about:
            claim["about"] = about
        ratings = [{"@type": "Rating", "reviewAspect": "verdict", "alternateName": c["verdict"] or "Unrated",
                    "ratingExplanation": c["reasons"] or None}]
        if c["certainty"] in CERTAINTY:
            ratings.append({"@type": "Rating", "reviewAspect": "certainty", "alternateName": c["certainty"],
                            "ratingValue": CERTAINTY[c["certainty"]], "worstRating": 1, "bestRating": 4})
        review = {"@id": f"#{cid}", "@type": "ClaimReview", "claimReviewed": c["statement"], "itemReviewed": claim,
                  "reviewRating": ratings, "reviewBody": c["evidence_text"] if "evidence_text" in have else None,
                  "keywords": [c["contrast"]] if c["contrast"] else None,
                  "isBasedOn": {"@id": "claim_store.sqlite"}}
        graph.append(_clean(review))
    con.close()
    return {"@context": {"@vocab": "https://schema.org/"}, "@graph": graph}


def _clean(x):
    if isinstance(x, dict):
        return {k: _clean(v) for k, v in x.items() if v not in (None, "", [], {})}
    if isinstance(x, list):
        return [_clean(v) for v in x]
    return x


def _authors():
    """Authors from CITATION.cff (family-names / given-names / orcid), without a YAML parser."""
    text = (ROOT / "CITATION.cff").read_text()
    block = text.split("\nauthors:", 1)[1] if "\nauthors:" in text else ""
    block = re.split(r"\n[A-Za-z][\w-]*:", block, maxsplit=1)[0]
    people = []
    for entry in re.split(r"\n\s*-\s+", "\n" + block)[1:]:
        f = dict(re.findall(r'([\w-]+):\s*"?([^"\n]+)"?', entry))
        name = " ".join(x for x in (f.get("given-names"), f.get("family-names")) if x)
        if name:
            people.append({"name": name, "orcid": f.get("orcid", "").strip()})
    return people


def ro_crate(out, info):
    out = Path(out)
    authors = _authors()
    people = []
    for i, a in enumerate(authors, 1):
        pid = a["orcid"] if a["orcid"].startswith("https://orcid.org/") else f"#author-{i}"
        people.append({"@id": pid, "@type": "Person", "name": a["name"]})
    files, dirs = [], []
    for p in sorted(out.rglob("*")):
        rel = p.relative_to(out).as_posix()
        if rel == "ro-crate-metadata.json":
            continue
        if p.is_dir():
            dirs.append({"@id": rel + "/", "@type": "Dataset", "name": rel})
        else:
            files.append({"@id": rel, "@type": "File", "name": p.name, "encodingFormat": FORMATS.get(p.suffix, "text/plain"),
                          "contentSize": str(p.stat().st_size)})
    files += [{"@id": f, "@type": "File", "name": f, "encodingFormat": FORMATS.get(Path(f).suffix, "text/plain")}
              for f in LATER if not (out / f).exists()]
    top = [{"@id": e["@id"]} for e in dirs + files if "/" not in e["@id"].rstrip("/")]
    for d in dirs:
        d["hasPart"] = [{"@id": e["@id"]} for e in dirs + files
                        if e["@id"].rstrip("/").rsplit("/", 1)[0] + "/" == d["@id"] and e is not d]
    root = {"@id": "./", "@type": "Dataset", "name": info.get("name", "Claim store release"),
            "description": info.get("description", ""), "version": info.get("version", ""),
            "datePublished": info.get("date", ""), "license": {"@id": CC_BY},
            "author": [{"@id": p["@id"]} for p in people], "keywords": info.get("keywords", []), "hasPart": top}
    graph = [{"@id": "ro-crate-metadata.json", "@type": "CreativeWork",
              "conformsTo": {"@id": "https://w3id.org/ro/crate/1.1"}, "about": {"@id": "./"}},
             _clean(root), {"@id": CC_BY, "@type": "CreativeWork", "name": "CC BY 4.0"}, *people, *dirs, *files]
    return {"@context": "https://w3id.org/ro/crate/1.1/context", "@graph": graph}


def write(out, info):
    out = Path(out)
    (out / "claims.jsonld").write_text(json.dumps(claims_jsonld(out / "claim_store.sqlite"), ensure_ascii=False) + "\n")
    (out / "ro-crate-metadata.json").write_text(json.dumps(ro_crate(out, info), indent=1, ensure_ascii=False) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("release_dir")
    a = ap.parse_args()
    out = Path(a.release_dir)
    info = {}
    if (out / "release.json").exists():
        m = json.loads((out / "release.json").read_text())
        info = {k: m.get(k, "") for k in ("name", "version", "date")}
    info["description"] = cfg.describe()["description"]
    write(out, info)
    print(f"{out}: claims.jsonld, ro-crate-metadata.json")


if __name__ == "__main__":
    main()

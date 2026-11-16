"""Linked-data exports of a release (scripts/59d_linked_data.py; COMMUNITY_PLAN step 2.5).

Standard library only:  python3 -I tests/test_linked_data.py   or   pytest tests/
"""
import importlib.util
import json
import shutil
import sqlite3
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixture" / "claim_store.sqlite"
SCHEMA_ORG = {"about", "alternateName", "appearance", "bestRating", "claimReviewed", "datePublished", "identifier",
              "inDefinedTermSet", "isBasedOn", "isPartOf", "itemReviewed", "keywords", "name", "ratingExplanation",
              "ratingValue", "reviewAspect", "reviewBody", "reviewRating", "termCode", "text", "worstRating"}


def _module():
    spec = importlib.util.spec_from_file_location("linked_data", ROOT / "scripts" / "59d_linked_data.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _terms(x, out):
    if isinstance(x, dict):
        for k, v in x.items():
            if not k.startswith("@"):
                out.add(k)
            _terms(v, out)
    elif isinstance(x, list):
        for v in x:
            _terms(v, out)
    return out


def test_claims_jsonld_covers_every_claim():
    doc = _module().claims_jsonld(FIXTURE)
    n = sqlite3.connect(f"file:{FIXTURE}?mode=ro", uri=True).execute("SELECT count(*) FROM claim").fetchone()[0]
    assert doc["@context"] == {"@vocab": "https://schema.org/"}
    assert len(doc["@graph"]) == n
    assert _terms(doc, set()) <= SCHEMA_ORG  # only the schema.org properties we checked by hand
    r = doc["@graph"][0]
    assert r["@type"] == "ClaimReview" and r["itemReviewed"]["@type"] == "Claim"
    assert r["itemReviewed"]["appearance"]["@id"].startswith("https://doi.org/")
    assert any(t.get("@id", "").startswith("http://purl.obolibrary.org/obo/") for t in r["itemReviewed"].get("about", []))


def test_ro_crate_lists_the_release(tmp_path=None):
    ld = _module()
    d = Path(tmp_path or tempfile.mkdtemp())
    shutil.copy(FIXTURE, d / "claim_store.sqlite")
    (d / "tables").mkdir()
    (d / "tables" / "claim.tsv").write_text("claim_id\n")
    ld.write(d, {"name": "test", "version": "0.0.1", "date": "2026-01-01"})
    crate = json.loads((d / "ro-crate-metadata.json").read_text())
    ids = {n["@id"]: n for n in crate["@graph"]}
    assert ids["ro-crate-metadata.json"]["about"] == {"@id": "./"}
    root = ids["./"]
    assert root["license"] == {"@id": ld.CC_BY} and root["author"]
    assert {"claim_store.sqlite", "claims.jsonld", "tables/", "release.json", "SHA256SUMS"} <= {p["@id"] for p in root["hasPart"]}
    assert ids["tables/"]["hasPart"] == [{"@id": "tables/claim.tsv"}]
    assert ids["claim_store.sqlite"]["encodingFormat"] == "application/vnd.sqlite3"


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print("PASS", name)

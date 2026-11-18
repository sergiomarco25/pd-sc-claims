"""Build a public release of the claim store (claims/COMMUNITY_PLAN.md, phase 4).

From the full store (config STORE_DB) it writes results/releases/v<VERSION>/:

  claim_store.sqlite          public store: private datasets, private papers' claims and findings computed on private
                              data removed (config [release]); findings re-derived without private data get their
                              public tests (59b_public_findings.py); every claim re-rated on what remains
  tables/<table>.tsv          every table as TSV
  knowledge_graph.graphml     graph of the public store (Cytoscape, Gephi, networkx)
  knowledge_graph.json        the same as node-link JSON
  LICENSE-DATA.md, CITATION.cff, RATING.md, schema/   licence, how to cite, rating rules, claim schemas
  CHANGES.md, changes.tsv     per-claim changes against the previous release (with --previous; 59c_release_diff.py)
  claims.jsonld               claims as schema.org ClaimReview JSON-LD (59d_linked_data.py)
  ro-crate-metadata.json      RO-Crate 1.1 description of the release folder (59d_linked_data.py)
  release.json                version, date, code commit, rating and schema versions, what was removed, file hashes
  SHA256SUMS                  checksums of every file

Re-derived findings' tests and ratings are recomputed with scripts/59b_public_findings.py and 52_rating.py (need
pandas and scipy: the analysis environment, or any Python given with --rating-python). If a re-derived finding cannot
be rebuilt it is removed; without a successful re-rating the release is marked "public": false and must not be
published. Uploading (Zenodo, Lamin) is a separate step.

Usage: python3 -I scripts/59_release.py VERSION [--previous OLD_STORE] [--out DIR] [--rating-python PATH] [--force]
"""
import argparse
import csv
import hashlib
import json
import re
import shutil
import sqlite3
import subprocess
import sys
import tempfile
from datetime import date, datetime, timezone
from pathlib import Path

from claimstore import audit, kg, linked_data, release_diff
from claimstore import config as cfg

ROOT = cfg.REPO

RATING_COLS = ["verdict", "certainty", "rating", "family", "start", "k_cohorts", "est", "lo", "hi", "i2", "concordance",
               "z_comb", "detect", "donors", "rob", "inconsistency", "indirectness", "imprecision", "dependence",
               "upgrade", "reasons", "not_assessed", "scope_tested", "scope_held", "caveats"]
OUTCOMES = ["reproduced", "direction_only", "not_reproduced", "contradicted"]


def _cols(con, table):
    return [r[1] for r in con.execute(f"PRAGMA table_info('{table}')")]


def _ph(n):
    return ",".join("?" * n) or "''"


def drop_claims(con, ids):
    """Delete claims and everything keyed on them (tests, relations, signature links)."""
    ids = sorted(ids)
    con.execute("CREATE TEMP TABLE IF NOT EXISTS dropc(claim_id TEXT PRIMARY KEY)")
    con.execute("DELETE FROM dropc")
    con.executemany("INSERT INTO dropc VALUES (?)", [(c,) for c in ids])
    for (t,) in con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'").fetchall():
        if "claim_id" in _cols(con, t):
            con.execute(f"DELETE FROM '{t}' WHERE claim_id IN (SELECT claim_id FROM dropc)")
    con.execute("DELETE FROM claim_relation WHERE claim_a IN (SELECT claim_id FROM dropc) "
                "OR claim_b IN (SELECT claim_id FROM dropc)")
    con.commit()


def make_public(con, private_datasets, private_papers, separable_findings=(), rederived_findings=()):
    """Remove private material in place; return what was removed. Findings with a test in a private dataset are
    left out unless listed as separable (config [release] separable_findings: they keep their other datasets) or as
    re-derived (rederived_findings: their tests are then replaced by 59b_public_findings.py)."""
    pd_, pp = sorted(private_datasets), sorted(private_papers)
    removed = {}
    drop = {r[0] for r in con.execute(f"SELECT claim_id FROM claim WHERE paper IN ({_ph(len(pp))})", pp)}
    removed["claims of private papers"] = len(drop)
    fin = {r[0] for r in con.execute(
        f"SELECT DISTINCT claim_id FROM claim_dataset WHERE claim_id LIKE 'F-%' AND dataset IN ({_ph(len(pd_))})", pd_)}
    kept, redo = fin & set(separable_findings), fin & set(rederived_findings)
    fin -= kept | redo
    removed["findings depending on private data"] = len(fin)
    removed["findings kept without their private-dataset tests"] = len(kept)
    removed["findings re-derived without private data"] = len(redo)
    drop |= fin
    removed["tests in private datasets"] = con.execute(
        f"SELECT count(*) FROM claim_dataset WHERE dataset IN ({_ph(len(pd_))})", pd_).fetchone()[0]
    drop_claims(con, drop)
    tables = [r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")]
    for t in tables:
        if "dataset" in _cols(con, t) and pd_:
            con.execute(f"DELETE FROM '{t}' WHERE dataset IN ({_ph(len(pd_))})", pd_)
    if pp:
        con.execute(f"DELETE FROM signature WHERE paper IN ({_ph(len(pp))})", pp)
        removed["private papers in the literature inventory"] = con.execute(
            f"DELETE FROM source WHERE paper_id IN ({_ph(len(pp))})", pp).rowcount
    con.execute("DELETE FROM claim_signature WHERE signature_id NOT IN (SELECT signature_id FROM signature)")
    con.execute("DELETE FROM signature WHERE signature_id NOT IN (SELECT signature_id FROM claim_signature)")
    con.execute("DELETE FROM signature_gene WHERE signature_id NOT IN (SELECT signature_id FROM signature)")
    # Our evidence notes on consensus statements may describe private results; questions come from the private atlas
    # (public questions go to the GitHub issue form). The legacy tier summary counted private tests: superseded by
    # the rating, so it is dropped rather than recomputed.
    con.execute("UPDATE consensus_statement SET our_evidence = ''")
    removed["questions"] = con.execute("SELECT count(*) FROM question").fetchone()[0]
    con.execute("DELETE FROM question")
    con.execute("DROP TABLE IF EXISTS claim_summary")
    con.execute("UPDATE claim SET tier = NULL, rank = NULL")
    con.commit()
    return removed


def redact(con, private_ids):
    """After re-rating: clear the source dataset of claims made on private data, and drop every sentence of a claim's
    text fields that names a private dataset or paper (curation notes can quote private results). Returns counts."""
    ids = sorted(private_ids)
    if not ids:
        return {}
    pat = re.compile(r"\b(" + "|".join(map(re.escape, ids)) + r")\b")
    n_src = con.execute(f"UPDATE claim SET source_dataset = NULL WHERE source_dataset IN ({_ph(len(ids))})", ids).rowcount
    keep = {"claim_id", "paper", "source_dataset"}
    cols = [c for c in _cols(con, "claim") if c not in keep]
    n_txt = 0
    for row in con.execute(f"SELECT claim_id, {', '.join(f'[{c}]' for c in cols)} FROM claim").fetchall():
        for c, v in zip(cols, row[1:]):
            if isinstance(v, str) and pat.search(v):
                kept = [x for x in re.split(r"(?<=\.)\s+(?=[A-Z(])", v) if not pat.search(x)]
                con.execute(f"UPDATE claim SET [{c}] = ? WHERE claim_id = ?", (" ".join(kept) or None, row[0]))
                n_txt += 1
    con.commit()
    return {"claims with a private source dataset cleared": n_src, "claim text fields redacted": n_txt}


def private_mentions(db):
    """Every text cell of every table that still names a private dataset, paper or module (scripts/61_audit_public.py
    patterns), as {table.column: count}."""
    pats = audit.patterns()
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    found = {}
    for (t,) in con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'").fetchall():
        cols = _cols(con, t)
        for row in con.execute(f"SELECT * FROM '{t}'"):
            for c, v in zip(cols, row):
                if isinstance(v, str) and any(p.search(v) for _, p in pats):
                    found[f"{t}.{c}"] = found.get(f"{t}.{c}", 0) + 1
    con.close()
    return found


def _engine_cmd(python, module, *args):
    """Command running claimstore.<module>.main() with another Python (the analysis environment)."""
    code = (f"import sys; sys.path.insert(0, {str(ROOT / 'src')!r}); sys.argv[0] = {module!r}; "
            f"from claimstore.{module} import main; main()")
    return [str(python), "-I", "-c", code, *map(str, args)]


def rederive(public_db, python):
    """Run claimstore.public_findings on the public store; return its summary {finding id: ...}."""
    r = subprocess.run(_engine_cmd(python, "public_findings", "--db", public_db), capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError("re-derived findings not rebuilt:\n" + r.stderr[-3000:])
    return json.loads(r.stdout.strip().splitlines()[-1])


def rerate(public_db, python):
    """Run claimstore.rating on the public store; return {claim_id: row} and the rating version."""
    with tempfile.TemporaryDirectory() as tmp:
        r = subprocess.run(_engine_cmd(python, "rating", "--db", public_db, "--out", tmp),
                           capture_output=True, text=True)
        if r.returncode != 0:
            raise RuntimeError("re-rating failed:\n" + r.stderr[-3000:])
        with open(Path(tmp) / "claim_rating.tsv", newline="") as fh:
            rows = {row["claim_id"]: row for row in csv.DictReader(fh, delimiter="\t")}
    version = next(iter(rows.values())).get("rating_version", "") if rows else ""
    return rows, version


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def apply_ratings(con, rows):
    """Write the new ratings into the claim table and rebuild each claim's evidence text from public tests only."""
    have = set(_cols(con, "claim"))
    cols = [c for c in RATING_COLS if c in have]
    counts = {}
    for cid, ds, role, outcome in con.execute("SELECT claim_id, dataset, role, outcome FROM claim_dataset"):
        if role == "source":
            continue
        k = counts.setdefault(cid, {"datasets": set(), **{o: 0 for o in OUTCOMES}})
        k["datasets"].add(ds)
        if outcome in OUTCOMES:
            k[outcome] += 1
    for (cid,) in con.execute("SELECT claim_id FROM claim").fetchall():
        rt = rows.get(cid, {})
        vals = [(_num(rt.get(c)) if c not in ("verdict", "certainty", "rating", "family", "start", "reasons",
                                                "not_assessed", "scope_tested", "scope_held", "caveats")
                 else (rt.get(c) or None)) for c in cols]
        k = counts.get(cid)
        text = ""
        if k:
            text = (f"Other datasets (public release): tested in {len(k['datasets'])} ({k['reproduced']} reproduced, "
                    f"{k['direction_only']} direction only, {k['not_reproduced']} not reproduced, "
                    f"{k['contradicted']} contradicted).")
        if rt:
            est = _num(rt.get("est"))
            eff = (f" Pooled effect (approx. SMD, positive = as claimed) {est:.2f}, 95% CI {_num(rt['lo']):.2f} to "
                   f"{_num(rt['hi']):.2f}." if est is not None else "")
            text += (f" Rating: {rt['verdict']}, {rt['certainty'].lower()} certainty "
                     f"({rt.get('k_cohorts') or 0} independent cohorts).{eff} Concerns: {rt.get('reasons', '')}.")
        con.execute(f"UPDATE claim SET {', '.join(c + ' = ?' for c in cols)}, evidence_text = ? WHERE claim_id = ?",
                    vals + [text.strip(), cid])
    con.commit()


def export_tables(db, out):
    out.mkdir(parents=True, exist_ok=True)
    con = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    for (t,) in con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'"):
        cur = con.execute(f"SELECT * FROM '{t}'")
        with open(out / f"{t}.tsv", "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow([d[0] for d in cur.description])
            w.writerows(cur)
    con.close()


def git_commit():
    try:
        r = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True, timeout=30)
        dirty = subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain"], capture_output=True, text=True,
                               timeout=60).stdout.strip()
        return r.stdout.strip() + ("-dirty" if dirty else "") if r.returncode == 0 else "unknown"
    except (OSError, subprocess.SubprocessError):
        return "unknown"


def sha256(p):
    h = hashlib.sha256()
    with open(p, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("version", help="semantic version, e.g. 1.0.0")
    ap.add_argument("--out", help="output folder (default: results/releases/v<VERSION>)")
    ap.add_argument("--db", default=str(cfg.STORE_DB), help="full store (default from config/project.toml)")
    ap.add_argument("--rating-python", default=str(cfg.PYTHON), help="Python with pandas and scipy for 52_rating.py")
    ap.add_argument("--force", action="store_true", help="overwrite an existing release folder")
    ap.add_argument("--previous", help="store of the previous release, for the per-claim changelog (CHANGES.md)")
    a = ap.parse_args()

    out = Path(a.out) if a.out else cfg.RESULTS / "releases" / f"v{a.version}"
    if out.exists():
        if not a.force:
            sys.exit(f"{out} exists (use --force to overwrite)")
        shutil.rmtree(out)
    out.mkdir(parents=True)
    db = out / "claim_store.sqlite"
    src = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
    con = sqlite3.connect(db)
    src.backup(con)
    src.close()
    removed = make_public(con, cfg.PRIVATE_DATASETS, cfg.PRIVATE_PAPERS, cfg.SEPARABLE_FINDINGS,
                          cfg.REDERIVED_FINDINGS)

    public, rating_version, rating_error = True, "", ""
    try:
        rederived = rederive(db, a.rating_python)
    except Exception as e:  # noqa: BLE001  (their old tests used private data: remove them, mark not public)
        redo = sorted(cfg.REDERIVED_FINDINGS)
        drop_claims(con, redo)
        rederived, public, rating_error = {"failed, removed": redo}, False, str(e)[-1500:]
        print(f"WARNING: re-derived findings removed, release marked not public:\n{rating_error}", file=sys.stderr)
    try:
        rows, rating_version = rerate(db, a.rating_python)
        apply_ratings(con, rows)
    except Exception as e:  # noqa: BLE001  (the release is then marked not public)
        public, rating_error = False, (rating_error + "\n" + str(e)[-1500:]).strip()
        print(f"WARNING: ratings not recomputed, release marked not public:\n{rating_error}", file=sys.stderr)
    removed.update(redact(con, cfg.PRIVATE_DATASETS | cfg.PRIVATE_PAPERS))
    built = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    con.execute("UPDATE build_info SET [commit] = ?, built = ?, n_claims = (SELECT count(*) FROM claim), "
                "n_datasets = (SELECT count(*) FROM dataset)", (f"release v{a.version}", built))
    con.commit()
    con.execute("VACUUM")
    counts = {t: con.execute(f"SELECT count(*) FROM '{t}'").fetchone()[0]
              for t in ("claim", "claim_dataset", "dataset", "consensus_statement", "background")}
    verdicts = dict(con.execute("SELECT verdict, count(*) FROM claim GROUP BY verdict ORDER BY 2 DESC").fetchall())
    con.close()
    leftover = private_mentions(db)
    if leftover:
        public = False
        rating_error = (rating_error + f"\nprivate names left in: {leftover}").strip()
        print(f"WARNING: private names left in {leftover}; release marked not public", file=sys.stderr)

    export_tables(db, out / "tables")
    g = kg.build(db)
    g.write_graphml(out / "knowledge_graph.graphml")
    (out / "knowledge_graph.json").write_text(json.dumps(g.to_json(), default=str))
    for f in ("LICENSE-DATA.md", "CITATION.cff"):
        shutil.copy2(ROOT / f, out / f)
    changes = None
    if a.previous:
        rows, labels = release_diff.diff(a.previous, db)
        changes = {"previous": labels[0], "counts": release_diff.write(rows, labels, out)}
    about = cfg.describe()
    linked_data.write(out, {"name": about["name"], "version": a.version, "date": str(date.today()),
                            "description": about["description"], "keywords": about["keywords"]})
    shutil.copy2(ROOT / "claims" / "RATING.md", out / "RATING.md")
    shutil.copytree(ROOT / "claims" / "schema", out / "schema", ignore=shutil.ignore_patterns("*.py", "__pycache__"))

    files = sorted(p for p in out.rglob("*") if p.is_file())
    schema_v = json.loads((ROOT / "claims" / "schema" / "claim.schema.json").read_text()).get("version", "")
    manifest = {
        "name": about["name"], "version": a.version, "date": str(date.today()), "built": built,
        "public": public, "code_commit": git_commit(), "rating_version": rating_version, "rating_error": rating_error,
        "claim_schema_version": schema_v, "licence": cfg.PROJECT.get("licence", ""),
        # counts, not names: a public release does not name what it leaves out (names stay in config/project.toml)
        "excluded": {"private_datasets": len(cfg.PRIVATE_DATASETS), "private_papers": len(cfg.PRIVATE_PAPERS),
                     "separable_findings": sorted(cfg.SEPARABLE_FINDINGS),
                     "rederived_findings": rederived, "removed": removed},
        "counts": counts, "verdicts": verdicts, "changes": changes, "private_mentions": leftover,
        "files": {str(p.relative_to(out)): {"bytes": p.stat().st_size, "sha256": sha256(p)} for p in files},
    }
    (out / "release.json").write_text(json.dumps(manifest, indent=2) + "\n")
    with open(out / "SHA256SUMS", "w") as fh:
        for p in sorted(p for p in out.rglob("*") if p.is_file() and p.name != "SHA256SUMS"):
            fh.write(f"{sha256(p)}  {p.relative_to(out)}\n")
    print(f"{out}: {counts['claim']} claims, {counts['claim_dataset']} tests, {counts['dataset']} datasets; "
          f"removed {removed}; public={public}")
    if not public:
        sys.exit(1)


if __name__ == "__main__":
    main()

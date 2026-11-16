"""Export the public repository (config [release] public_repository) from this private one.

Usage: python3 -I scripts/62_public_repo.py OUT --release DIR [--atlas FILE] [--force] [--allow-dirty]
Writes a fresh tree at OUT (no history) from the files committed in this repository (the working tree must be clean):
  left out   manuscript, reports, figures, internal plans and history, local tool settings, outputs of the full store
             (claims/store, merged claims, evidence tables, the question log), every file named after a private dataset
             or paper, and analyses that read private data (EXCLUDE)
  filtered   tables and YAML listing claims or datasets: rows and list items of private datasets and papers removed;
             findings that are not in the public release removed; the question log emptied
  redacted   curation notes, claim files and docs: sentences naming a private dataset or paper are dropped, as the
             release does (59_release.redact)
  replaced   README.md and CHANGELOG.md (public versions below), CITATION.cff (the public repository's URL),
             config/project.toml (private datasets out of the pipeline's dataset list)
  added      release/ (a public release from 59_release.py: store, tables, graph, linked data, checksums) and
             docs/index.html (the public atlas, for GitHub Pages)
Every text file is then scanned. A private name outside the code and config files (CODE) fails the export; those code
mentions, and wider words that can point at private material (REVIEW), are written to OUT/../<name>_review.tsv for a
person to read before the repository is made public. Needs pyyaml.
"""
import argparse
import csv
import fnmatch
import gzip
import io
import json
import re
import shutil
import subprocess
import sys
from datetime import date
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep import config as cfg  # noqa: E402

PRIV = sorted(cfg.PRIVATE_DATASETS | cfg.PRIVATE_PAPERS)
PAT = re.compile(r"\b(" + "|".join(map(re.escape, PRIV)) + r")\b") if PRIV else re.compile(r"(?!x)x")
MODULES = [re.compile(re.escape(m), re.I) for m in sorted(cfg.PRIVATE_MODULES)]
REVIEW = re.compile(r"GLIAD|\bMGH\b|Barba-Reyes(?! 2025)|revision/|unpublished|preprint", re.I)

EXCLUDE = [
    "manuscript/*", "reports/*", "figures/*", "PLAN.md", "CHANGELOG.md", "README.md", ".mcp.json", ".claude/*",
    "claims/store/*", "claims/claims.yaml", "claims/claims.tsv", "claims/claim_evidence.tsv",
    "claims/consensus_summary.tsv", "claims/source_discrepancies.tsv", "claims/COMMUNITY_PLAN.md", "claims/ROADMAP.md",
    "scripts/run_phase*.sh", "scripts/subgroups/0[4-9]_*", "scripts/survival/*", "scripts/38_spatial_stress.py",
    # findings entry (holds the private findings' statements), internal R1/R2 reports, manuscript figures and documents
    "scripts/53_findings.py", "scripts/r1_template.html", "scripts/r2_template.html", "scripts/02_build_r1_report.py",
    "scripts/06_build_r2_report.py", "scripts/18_figures.py", "scripts/19_supplementary_tables.py",
    "scripts/20_references.py", "scripts/21_build_docx.mjs", "scripts/31_story_figures.py", "scripts/32_build_pdf.py",
    "scripts/36_claim_figure.py", "scripts/40_map_figures.py", "scripts/49_video.py",
]
# code and config: a private dataset id here is a reference (a list entry, a label), not data; listed for review
CODE = ["*.py", "*.R", "*.sh", "*.sbatch", "*.toml", "*.yml", "*.html", "*.js", "workflow/Snakefile", ".gitignore"]
TEXT = {".md", ".tsv", ".txt", ".yaml", ".yml", ".py", ".R", ".sh", ".toml", ".html", ".js", ".cff", ".json", ".sbatch",
        ".csv", ".jsonl", ""}


def git(*args):
    return subprocess.run(["git", "-C", str(ROOT), *args], capture_output=True, text=True, check=True).stdout


def matches(path, globs):
    return any(fnmatch.fnmatch(path, g) for g in globs)


def named_private(path):
    return any(p in path for p in PRIV)


def redact_text(s):
    """Drop sentences that name a private dataset or paper (sentence: up to a full stop before a capital or '(')."""
    if not PAT.search(s):
        return s
    lines = []
    for line in s.split("\n"):
        if not PAT.search(line):
            lines.append(line)
            continue
        indent = re.match(r"\s*(?:[-*]\s+)?", line).group(0)
        kept = [x for x in re.split(r"(?<=\.)\s+(?=[A-Z(`])", line[len(indent):]) if not PAT.search(x)]
        if kept:
            lines.append(indent + " ".join(kept))
    return "\n".join(lines)


def redact_obj(x):
    if isinstance(x, str):
        return redact_text(x)
    if isinstance(x, list):
        return [redact_obj(v) for v in x if not (isinstance(v, str) and PAT.fullmatch(v))]
    if isinstance(x, dict):
        return {k: redact_obj(v) for k, v in x.items()}
    return x


def private_row(cells):
    key = cells[0] if cells else ""
    return bool(PAT.fullmatch(key) or any(key.startswith(p + "-") for p in PRIV))


def filter_tsv(text):
    """Rows keyed by a private dataset, paper or claim are dropped; ';'-lists lose private items; other cells naming
    one are redacted. Returns (text, rows dropped)."""
    rows = list(csv.reader(io.StringIO(text), delimiter="\t"))
    out, dropped = rows[:1], 0
    for r in rows[1:]:
        if private_row(r):
            dropped += 1
            continue
        out.append([";".join(x for x in c.split(";") if not PAT.search(x)) if ";" in c else redact_text(c) for c in r])
    buf = io.StringIO()
    csv.writer(buf, delimiter="\t", lineterminator="\n").writerows(out)
    return buf.getvalue(), dropped


def dump_yaml(obj):
    return yaml.safe_dump(obj, sort_keys=False, allow_unicode=True, width=120)


def export_file(rel, src, dst, public_findings, log):
    data = src.read_bytes()
    if rel == "claims/questions.yaml":  # the private atlas's question log; public questions come as GitHub issues
        data = b"[]\n"
        log.append((rel, "emptied"))
    elif rel == "claims/findings.yaml":
        items = yaml.safe_load(data) or []
        keep = [f for f in items if f["id"] in public_findings]
        for f in keep:
            flt = f.get("evidence", {}).get("filter") or {}
            for k, v in flt.items():
                if isinstance(v, list):
                    flt[k] = [x for x in v if not (isinstance(x, str) and PAT.fullmatch(x))]
        log.append((rel, f"findings kept {len(keep)} of {len(items)}"))
        data = dump_yaml(redact_obj(keep)).encode()
    elif rel == "config/project.toml":
        s = data.decode()
        for p in sorted(cfg.PRIVATE_DATASETS):
            s = re.sub(rf'"{re.escape(p)}",\s*', "", s)  # pipeline dataset list
            s = re.sub(rf"(?m)^{re.escape(p)}\s*=.*\n", "", s)  # its memory line
        log.append((rel, "private datasets removed from the pipeline list"))
        data = s.encode()
    elif rel == "CITATION.cff":
        s = re.sub(r'(?m)^repository-code:.*$', f'repository-code: "{cfg.PUBLIC_REPOSITORY}"', data.decode())
        data = s.encode()
    elif rel.endswith(".tsv.gz"):
        text = gzip.decompress(data).decode()
        if PAT.search(text):
            rows = list(csv.reader(io.StringIO(text), delimiter="\t"))
            hdr, n0 = rows[0], len(rows) - 1
            idx = [i for i, h in enumerate(hdr) if h in ("paper", "claim_ids", "signature_id")]
            rows = rows[:1] + [r for r in rows[1:] if not any(PAT.search(r[i]) for i in idx if i < len(r))]
            buf = io.StringIO()
            csv.writer(buf, delimiter="\t", lineterminator="\n").writerows(rows)
            data = gzip.compress(buf.getvalue().encode(), mtime=0)
            log.append((rel, f"rows of private papers dropped: {n0 - (len(rows) - 1)}"))
    elif rel.endswith(".tsv") and PAT.search(data.decode(errors="ignore")):
        text, n = filter_tsv(data.decode())
        log.append((rel, f"private rows dropped: {n}"))
        data = text.encode()
    elif rel.endswith((".yaml", ".yml")) and not matches(rel, CODE) and PAT.search(data.decode(errors="ignore")):
        obj = yaml.safe_load(data)
        if isinstance(obj, list):
            n0 = len(obj)
            obj = [e for e in obj if not (isinstance(e, dict) and any(isinstance(e.get(k), str) and PAT.fullmatch(e[k])
                                                                       for k in ("dataset", "paper", "id", "paper_id")))]
            log.append((rel, f"private entries dropped: {n0 - len(obj)}"))
        data = dump_yaml(redact_obj(obj)).encode()
        log.append((rel, "redacted"))
    elif Path(rel).suffix in {".md", ".txt"} and PAT.search(data.decode(errors="ignore")):
        data = redact_text(data.decode()).encode()
        log.append((rel, "redacted"))
    dst.parent.mkdir(parents=True, exist_ok=True)
    dst.write_bytes(data)
    shutil.copymode(src, dst)


README = """# PD single-cell claims: a re-tested evidence base

Published conclusions of human post-mortem brain Parkinson's disease single-cell and single-nucleus studies, extracted
as structured claims, each re-tested donor-level in its own data and in every other harmonised dataset, and rated for
certainty (verdict and GRADE-style certainty; rules in `claims/RATING.md`). The same store holds this study's own
findings, rated the same way, consensus statements and reference knowledge, and is served as an atlas, a knowledge
graph and a server for AI tools (Model Context Protocol).

**Status:** {status}

## What is here
| Path | What |
|---|---|
| `release/` | the current public release: `claim_store.sqlite`, every table as TSV, the knowledge graph (GraphML, JSON), claims as schema.org JSON-LD, an RO-Crate description, `release.json` and `SHA256SUMS` |
| `docs/index.html` | the atlas, a single static page built from the release (also served with GitHub Pages) |
| `claims/` | the claims of each published paper (`<paper>.yaml`, schema in `claims/schema/`), signatures, contexts, consensus statements, reference knowledge and the rating rules |
| `src/pdrep/`, `scripts/`, `workflow/` | the pipeline: harmonising datasets, donor-level tests, meta-analysis, rating, the store, the atlas and the release (Snakemake; `config/project.toml` holds every path) |
| `deploy/mcp/` | the MCP server: run it locally against `release/claim_store.sqlite` |
| `tests/` | tests on a small public test store (run by CI) |

One harmonised cohort whose own paper is not yet published is not part of this release: its data, its results and its
paper's claims are left out, and every rating is computed without it. It will be added after publication.

## Use
- Browse: open `docs/index.html`.
- Query: `sqlite3 release/claim_store.sqlite`, or the TSV files in `release/tables/`.
- AI tools: see `deploy/mcp/` to serve the store over MCP.
- Rebuild: the pipeline needs the raw datasets (accessions in `inventory/`); `snakemake -n` shows the plan.
- Ask a question of the evidence base: open an issue with the "Ask the evidence base" form.

## Licence and citation
Code: Apache-2.0 (`LICENSE`). Claims, store, tables and other data and text created here: CC BY 4.0 (`LICENSE-DATA.md`).
Harmonised data keep their sources' licences. Cite with `CITATION.cff`.
"""

CHANGELOG = """# Changelog

## {day} — First public export
- Exported from the private development repository with `scripts/62_public_repo.py` (code, published papers' claims,
  schema, rating rules and docs; private material left out or redacted). Release: {version}.
"""


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out")
    ap.add_argument("--release", required=True, help="a public release folder (59_release.py output, public: true)")
    ap.add_argument("--atlas", help="the public atlas page (51c_store_browser.py --public)")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--allow-dirty", action="store_true", help="preview only: export uncommitted changes too")
    a = ap.parse_args()
    out, rel_dir = Path(a.out), Path(a.release)
    if git("status", "--porcelain").strip() and not a.allow_dirty:
        sys.exit("the working tree has uncommitted changes: export from a clean checkout (or preview with --allow-dirty)")
    manifest = json.loads((rel_dir / "release.json").read_text())
    if not manifest.get("public"):
        sys.exit(f"{rel_dir} is not a public release")
    if out.exists():
        if not a.force:
            sys.exit(f"{out} exists (use --force)")
        shutil.rmtree(out)
    out.mkdir(parents=True)
    import sqlite3
    con = sqlite3.connect(f"file:{rel_dir / 'claim_store.sqlite'}?mode=ro", uri=True)
    public_findings = {r[0] for r in con.execute("SELECT claim_id FROM claim WHERE claim_id LIKE 'F-%'")}
    con.close()

    log, skipped = [], []
    files = git("ls-files").splitlines()
    if a.allow_dirty:  # new files not yet committed
        files += git("ls-files", "--others", "--exclude-standard").splitlines()
    for rel in files:
        if matches(rel, EXCLUDE) or named_private(rel):
            skipped.append(rel)
            continue
        export_file(rel, ROOT / rel, out / rel, public_findings, log)
    shutil.copytree(rel_dir, out / "release")
    if a.atlas:
        (out / "docs").mkdir(exist_ok=True)
        shutil.copy2(a.atlas, out / "docs" / "index.html")
    status = (f"release {manifest['version']} ({manifest['counts']['claim']} claims, {manifest['counts']['dataset']} "
              f"datasets, built {manifest['date']}).")
    (out / "README.md").write_text(README.format(status=status))
    (out / "CHANGELOG.md").write_text(CHANGELOG.format(day=date.today(), version=manifest["version"]))

    # audit
    fails, review = [], []
    for p in sorted(x for x in out.rglob("*") if x.is_file()):
        rel = p.relative_to(out).as_posix()
        if p.suffix not in TEXT or rel.startswith("release/") and p.suffix == ".sqlite":
            continue
        try:
            text = p.read_text()
        except UnicodeDecodeError:
            continue
        for i, line in enumerate(text.splitlines(), 1):
            priv = PAT.search(line) or any(m.search(line) for m in MODULES)
            if priv and not matches(rel, CODE):
                fails.append((rel, i, line.strip()[:160]))
            elif priv or REVIEW.search(line):
                review.append((rel, i, "private name in code" if priv else "review word", line.strip()[:160]))
    rep = out.parent / f"{out.name}_review.tsv"
    with open(rep, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(["file", "line", "kind", "text"])
        w.writerows(review)
    n_files = sum(1 for x in out.rglob("*") if x.is_file())
    size = sum(x.stat().st_size for x in out.rglob("*") if x.is_file()) / 1e6
    print(f"{out}: {n_files} files, {size:.1f} MB; left out {len(skipped)} files; changed {len(log)}")
    for rel, what in log:
        print(f"  {rel}: {what}")
    print(f"review list: {rep} ({len(review)} lines in {len({r[0] for r in review})} files)")
    if fails:
        print(f"FAILED: {len(fails)} private mentions outside code and config:")
        for f in fails[:40]:
            print("  ", *f)
        sys.exit(1)
    print("no private name outside code and config")


if __name__ == "__main__":
    main()

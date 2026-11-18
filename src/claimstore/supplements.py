"""Supplementary files of the claim-source papers: list them, or download them into raw/supplements/<paper>/.

Usage: python3 -I src/claimstore/cli.py supplements [--list] [--only REGEX] [paper_id ...]

PMCIDs come from inventory/papers.tsv: a `pmcid` column if the project has one, else the paper's `pmid` converted by
the NCBI ID converter. The files are those named as supplementary material in the paper's PMC full-text XML (Europe
PMC, or NCBI efetch when Europe PMC fails). Each is taken from the Europe PMC supplementary zip (open-access papers),
else from PMC's own copy of the article files, else from the Springer static-content server (DOIs 10.1038 / 10.1186 /
10.1007). Zip members are written by basename only.

--only keeps the listed files whose name matches REGEX in full (the ones approved for download).
--list downloads nothing: it prints each listed file with its size where a server reports one, so the downloads can be
approved first. Without it, writes inventory/supplements_manifest.tsv (paper, file, source, bytes, status). From the
reference project's scripts/07_fetch_supplements.py, which keeps its own PMCID table.
"""
import io
import json
import re
import sys
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path

import pandas as pd

from claimstore import config as cfg

ROOT = cfg.REPO
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
KEEP = re.compile(r"\.(xlsx?|csv|tsv|txt|zip|docx?|pdf|gz)$", re.I)
UA = {"User-Agent": "Mozilla/5.0 (research data retrieval)"}


def get(url, timeout=600):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
        return r.read()


def get_file(url, timeout=600):
    """A supplementary file's bytes; a web page in its place (an error or a bot check) is refused, never saved."""
    b = get(url, timeout)
    if b[:512].lstrip().lower().startswith((b"<!doctype html", b"<html")):
        raise ValueError("server sent a web page instead of the file")
    return b


def size(url, timeout=60):
    """Bytes a server reports for `url` (HEAD), or None when it does not say or does not have the file."""
    try:
        req = urllib.request.Request(url, headers=UA, method="HEAD")
        with urllib.request.urlopen(req, timeout=timeout) as r:
            n = r.headers.get("Content-Length")
            html = "html" in str(r.headers.get("Content-Type", "")).lower()
            return int(n) if n and int(n) > 0 and not html else None  # 0 or a web page: not the file
    except Exception:
        return None


def pmcids(papers):
    """paper_id -> PMCID ('' when the paper is not in PMC)."""
    inv = pd.read_csv(ROOT / "inventory" / "papers.tsv", sep="\t", dtype=str).fillna("").set_index("paper_id")
    out, todo = {}, {}
    for p in papers:
        row = inv.loc[p] if p in inv.index else pd.Series(dtype=str)
        if str(row.get("pmcid", "")).startswith("PMC"):
            out[p] = row["pmcid"]
        elif str(row.get("pmid", "")).strip().isdigit():
            todo[p] = row["pmid"].strip()
        else:
            out[p] = ""
    if todo:
        url = ("https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?format=json&tool=claimstore&ids="
               + ",".join(todo.values()))
        recs = {str(r.get("requested-id", "")): r.get("pmcid", "") or ""
                for r in json.loads(get(url, 60)).get("records", [])}
        out.update({p: recs.get(m, "") for p, m in todo.items()})
    return out


def listed(xml):
    """Supplementary file names and the DOI in a JATS full-text XML."""
    doi = re.search(r'<article-id pub-id-type="doi">([^<]+)<', xml)
    names = []
    for block in re.findall(r"<supplementary-material.*?</supplementary-material>", xml, re.S):
        names += re.findall(r'xlink:href="([^"]+)"', block)
    return sorted({n for n in names if KEEP.search(n)}), doi.group(1) if doi else ""


def full_text(pmcid):
    try:
        return get(f"{EPMC}/{pmcid}/fullTextXML").decode("utf-8", "replace")
    except Exception:  # Europe PMC answers 404 outside the open-access subset, sometimes 500; NCBI serves the same XML
        return get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=" + pmcid[3:]).decode(
            "utf-8", "replace")


def file_urls(pmcid, doi, name):
    """Where a single listed file can be fetched, in order of preference."""
    urls = [f"https://pmc.ncbi.nlm.nih.gov/articles/{pmcid}/bin/{urllib.parse.quote(name)}"]
    if doi.startswith(("10.1038", "10.1186", "10.1007")):
        urls.append("https://static-content.springer.com/esm/art%3A" + urllib.parse.quote(doi, safe="")
                    + "/MediaObjects/" + urllib.parse.quote(name))
    return urls


def _keep(names, only):
    return [n for n in names if only is None or re.fullmatch(only, n)]


def list_only(ids, only=None):
    rows = []
    for p, pmcid in ids.items():
        if not pmcid:
            rows.append({"paper": p, "pmcid": "", "file": "", "bytes": None, "where": "not in PMC"}); continue
        try:
            names, doi = listed(full_text(pmcid))
            names = _keep(names, only)
        except Exception as e:
            rows.append({"paper": p, "pmcid": pmcid, "file": "", "bytes": None, "where": f"XML failed: {e}"}); continue
        if not names:
            rows.append({"paper": p, "pmcid": pmcid, "file": "", "bytes": None, "where": "no supplementary files listed"})
        for n in names:
            b, where = None, ""
            for u in file_urls(pmcid, doi, n):
                b = size(u)
                if b is not None:
                    where = urllib.parse.urlparse(u).netloc
                    break
            rows.append({"paper": p, "pmcid": pmcid, "file": n, "bytes": b, "where": where or "size unknown"})
        time.sleep(0.5)
    t = pd.DataFrame(rows)
    with pd.option_context("display.max_rows", None, "display.width", 200, "display.max_colwidth", 70):
        print(t.to_string(index=False))
    known = t.bytes.dropna()
    print(f"\n{(t.file != '').sum()} files in {t.paper.nunique()} papers; {len(known)} with a size, "
          f"{known.sum() / 1e6:.1f} MB in total")
    return t


def fetch(ids, only=None):
    dest = cfg.RAW / "supplements"
    man_path = ROOT / "inventory" / "supplements_manifest.tsv"
    old = pd.read_csv(man_path, sep="\t") if man_path.exists() else pd.DataFrame()
    rows = []
    for p, pmcid in ids.items():
        if not pmcid:
            rows.append({"paper": p, "file": "", "source": "", "bytes": 0, "status": "not in PMC"}); continue
        out = dest / p
        out.mkdir(parents=True, exist_ok=True)
        try:
            names, doi = listed(full_text(pmcid))
            names = _keep(names, only)
        except Exception as e:
            rows.append({"paper": p, "file": "", "source": "", "bytes": 0, "status": f"XML failed: {e}"}); continue
        got = {f.name for f in out.iterdir()}
        present = set(names) & got
        if set(names) - got:
            try:
                z = zipfile.ZipFile(io.BytesIO(get(f"{EPMC}/{pmcid}/supplementaryFiles")))
                for m in z.namelist():
                    base = Path(m).name
                    if base in names and base not in got:
                        (out / base).write_bytes(z.read(m)); got.add(base)
                        rows.append({"paper": p, "file": base, "source": "europepmc", "bytes": (out / base).stat().st_size,
                                     "status": "ok"})
            except Exception as e:
                print(p, "europepmc zip failed:", e)
        for n in sorted(set(names) - got):
            for u in file_urls(pmcid, doi, n):
                try:
                    (out / n).write_bytes(get_file(u)); got.add(n)
                    rows.append({"paper": p, "file": n, "source": urllib.parse.urlparse(u).netloc,
                                 "bytes": (out / n).stat().st_size, "status": "ok"})
                    break
                except Exception as e:
                    print(p, n, urllib.parse.urlparse(u).netloc, "failed:", e)
            else:
                rows.append({"paper": p, "file": n, "source": "", "bytes": 0, "status": "missing"})
        for n in sorted(present):
            rows.append({"paper": p, "file": n, "source": "already present", "bytes": (out / n).stat().st_size,
                         "status": "ok"})
        print(p, f"{len(got & set(names))}/{len(names)} files", flush=True)
        time.sleep(1)
    new = pd.DataFrame(rows)
    man = pd.concat([old[~old.paper.isin(list(ids))] if len(old) else old, new]).sort_values(["paper", "file"])
    man.to_csv(man_path, sep="\t", index=False)
    print(man.groupby(["paper", "status"]).size().unstack(fill_value=0).to_string())
    return man


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    only_list = "--list" in argv
    only = argv[argv.index("--only") + 1] if "--only" in argv else None
    papers = [a for i, a in enumerate(argv) if a not in ("--list", "--only") and (i == 0 or argv[i - 1] != "--only")]
    if not papers:
        papers = pd.read_csv(ROOT / "inventory" / "papers.tsv", sep="\t", dtype=str).paper_id.tolist()
    ids = pmcids(papers)
    (list_only if only_list else fetch)(ids, only)
    return 0


if __name__ == "__main__":
    sys.exit(main())

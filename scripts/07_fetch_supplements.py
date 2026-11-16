"""Download the supplementary files of every claim-source paper.

Usage: python -I scripts/07_fetch_supplements.py [paper_id ...]

For each paper: list the supplementary files named in the Europe PMC full-text XML; download the Europe PMC
supplementary zip and keep the listed files; fetch any listed file still missing from the Springer static-content
server (DOIs 10.1038 / 10.1186 / 10.1007). Files go to raw/supplements/<paper>/ (one folder per paper; zip members
are written by basename only). Writes inventory/supplements_manifest.tsv (paper, file, source, bytes, status).
"""
from pathlib import Path
import io
import re
import sys
import time
import urllib.parse
import urllib.request
import zipfile
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

DEST = cfg.RAW / "supplements"
PMCID = {"adams2024": "PMC11361719", "barbareyes2025": "PMC12053221", "barbareyes2026": "PMC13142618",
         "chia2025": "PMC12359983", "clarence2025": "PMC12324657", "dehestani2024": "PMC11323592",
         "feleke2021": "PMC8357687", "jin2024": "PMC11389790", "kamath2022": "PMC9076534", "lee2023": "PMC10104466",
         "lee2026": "PMC13600946", "ma2025": "PMC12322057", "martirosyan2024": "PMC10799528", "menon2026": "PMC13001318",
         "prashant2024": "PMC11585549", "smajic2022": "PMC9050543", "wang2024": "PMC10780895", "xu2023": "PMC10046522",
         "zhu2024": "PMC12372474"}
EPMC = "https://www.ebi.ac.uk/europepmc/webservices/rest"
KEEP = re.compile(r"\.(xlsx?|csv|tsv|txt|zip|docx?|pdf|gz)$", re.I)
UA = {"User-Agent": "Mozilla/5.0 (research data retrieval)"}


def get(url, timeout=600):
    with urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=timeout) as r:
        return r.read()


def listed_files(pmcid):
    try:
        xml = get(f"{EPMC}/{pmcid}/fullTextXML").decode("utf-8", "replace")
    except Exception:  # Europe PMC sometimes answers 500; NCBI efetch serves the same JATS XML
        xml = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pmc&id=" + pmcid[3:]).decode("utf-8", "replace")
    doi = re.search(r'<article-id pub-id-type="doi">([^<]+)<', xml)
    names = []
    for block in re.findall(r"<supplementary-material.*?</supplementary-material>", xml, re.S):
        names += re.findall(r'xlink:href="([^"]+)"', block)
    return sorted({n for n in names if KEEP.search(n)}), doi.group(1) if doi else ""


def main(papers):
    papers = papers or list(PMCID)
    man_path = ROOT / "inventory" / "supplements_manifest.tsv"
    old = pd.read_csv(man_path, sep="\t") if man_path.exists() else pd.DataFrame()
    rows = []
    for p in papers:
        out = DEST / p
        out.mkdir(parents=True, exist_ok=True)
        try:
            names, doi = listed_files(PMCID[p])
        except Exception as e:
            rows.append({"paper": p, "file": "", "source": "", "bytes": 0, "status": f"XML failed: {e}"}); continue
        got = {f.name for f in out.iterdir()}
        if set(names) - got:
            try:
                z = zipfile.ZipFile(io.BytesIO(get(f"{EPMC}/{PMCID[p]}/supplementaryFiles")))
                for m in z.namelist():
                    base = Path(m).name
                    if base in names and base not in got:
                        (out / base).write_bytes(z.read(m)); got.add(base)
                        rows.append({"paper": p, "file": base, "source": "europepmc", "bytes": (out / base).stat().st_size, "status": "ok"})
            except Exception as e:
                print(p, "europepmc zip failed:", e)
        for n in sorted(set(names) - got):
            if doi.startswith(("10.1038", "10.1186", "10.1007")):
                url = ("https://static-content.springer.com/esm/art%3A" + urllib.parse.quote(doi, safe="")
                       + "/MediaObjects/" + n)
                try:
                    (out / n).write_bytes(get(url)); got.add(n)
                    rows.append({"paper": p, "file": n, "source": "springer", "bytes": (out / n).stat().st_size, "status": "ok"})
                    continue
                except Exception as e:
                    print(p, n, "springer failed:", e)
            rows.append({"paper": p, "file": n, "source": "", "bytes": 0, "status": "missing"})
        for n in sorted(set(names) & got - {r["file"] for r in rows if r["paper"] == p}):
            rows.append({"paper": p, "file": n, "source": "already present", "bytes": (out / n).stat().st_size, "status": "ok"})
        print(p, f"{len(got & set(names))}/{len(names)} files", flush=True)
        time.sleep(1)
    new = pd.DataFrame(rows)
    man = pd.concat([old[~old.paper.isin(papers)] if len(old) else old, new]).sort_values(["paper", "file"])
    man.to_csv(man_path, sep="\t", index=False)
    print(man.groupby(["paper", "status"]).size().unstack(fill_value=0).to_string())


if __name__ == "__main__":
    main(sys.argv[1:])

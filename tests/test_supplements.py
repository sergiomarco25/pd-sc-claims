"""Supplementary files (claimstore.supplements): parsing of the full-text XML and the download addresses; no network."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore import supplements as S  # noqa: E402

XML = """<article><front><article-meta><article-id pub-id-type="doi">10.1038/s41586-019-1195-2</article-id>
<supplementary-material id="S1"><media xlink:href="41586_2019_1195_MOESM4_ESM.xlsx"/></supplementary-material>
<supplementary-material id="S2"><media xlink:href="NIHMS1527283-supplement-2.pdf"/></supplementary-material>
<supplementary-material id="S3"><media xlink:href="movie.mp4"/></supplementary-material>
</article-meta></front></article>"""


def test_listed_keeps_tables_and_doi():
    names, doi = S.listed(XML)
    assert names == ["41586_2019_1195_MOESM4_ESM.xlsx", "NIHMS1527283-supplement-2.pdf"]
    assert doi == "10.1038/s41586-019-1195-2"


def test_file_urls_pmc_first_then_springer():
    urls = S.file_urls("PMC6865822", "10.1038/s41586-019-1195-2", "a b.xlsx")
    assert urls[0] == "https://pmc.ncbi.nlm.nih.gov/articles/PMC6865822/bin/a%20b.xlsx"
    assert urls[1].startswith("https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-019-1195-2/MediaObjects/")
    assert S.file_urls("PMC1", "10.1016/j.cell.2023.08.039", "x.xlsx") == [
        "https://pmc.ncbi.nlm.nih.gov/articles/PMC1/bin/x.xlsx"]


def test_web_page_is_not_saved_as_a_file(monkeypatch):
    monkeypatch.setattr(S, "get", lambda url, timeout=600: b"  <!DOCTYPE html><html>checking your browser</html>")
    try:
        S.get_file("https://example.org/x.xlsx")
    except ValueError:
        pass
    else:
        raise AssertionError("a web page was accepted as a file")
    monkeypatch.setattr(S, "get", lambda url, timeout=600: b"PK\x03\x04 table bytes")
    assert S.get_file("https://example.org/x.xlsx").startswith(b"PK")


def test_only_keeps_approved_files():
    names = ["NIHMS1-supplement-1.pdf", "NIHMS1-supplement-10.xls", "NIHMS1-supplement-3.xlsx", "x.sd01.xlsx"]
    assert S._keep(names, r"NIHMS1-supplement-(3|10)\.xlsx?") == ["NIHMS1-supplement-10.xls", "NIHMS1-supplement-3.xlsx"]
    assert S._keep(names, None) == names

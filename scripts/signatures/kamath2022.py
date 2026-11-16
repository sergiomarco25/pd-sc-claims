"""Gene signatures for kamath2022 claims (see claims/SIGNATURES.md).

Source: Supplementary Tables (41593_2022_1061_MOESM3_ESM.xlsx). The workbook is saved in
strict OOXML, which openpyxl cannot read, so a minimal XML reader is included below.
Only Supplementary Table 8 (MAST markers for the 10 DA subtypes, top 3,500 per subtype) holds
gene-level statistics; there is no PD-vs-control DE table and no glial-subtype marker table.

Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/kamath2022.py
"""
import re
import zipfile
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import yaml

PAPER = "kamath2022"
ROOT = Path(__file__).resolve().parents[2]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

SUPP = cfg.RAW / "supplements" / PAPER
XLSX = SUPP / "41593_2022_1061_MOESM3_ESM.xlsx"
OUT = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison",
        "gene", "rank", "stat", "padj", "source", "selection"]


# ---------------------------------------------------------------- strict-OOXML reader
def _strip(tag):
    return tag.split("}", 1)[1] if "}" in tag else tag


def _col_idx(ref):
    n = 0
    for ch in re.match(r"[A-Z]+", ref).group(0):
        n = n * 26 + (ord(ch) - 64)
    return n - 1


def read_strict_xlsx(path, sheets=None):
    """Return {sheet_name: DataFrame(header=None)} with cell values as strings/numbers."""
    z = zipfile.ZipFile(path)
    ss = []
    if "xl/sharedStrings.xml" in z.namelist():
        for si in ET.fromstring(z.read("xl/sharedStrings.xml")):
            ss.append("".join(t.text or "" for t in si.iter() if _strip(t.tag) == "t"))
    wb = ET.fromstring(z.read("xl/workbook.xml"))
    rels = {r.get("Id"): r.get("Target") for r in ET.fromstring(z.read("xl/_rels/workbook.xml.rels"))}
    out = {}
    for sh in wb.iter():
        if _strip(sh.tag) != "sheet":
            continue
        name = sh.get("name")
        if sheets and name not in sheets:
            continue
        rid = [v for k, v in sh.attrib.items() if _strip(k) == "id"][0]
        target = rels[rid].lstrip("/")
        target = target if target.startswith("xl/") else "xl/" + target
        rows, maxc = {}, 0
        for _, el in ET.iterparse(z.open(target)):
            if _strip(el.tag) != "c":
                continue
            ref, t, val = el.get("r"), el.get("t"), None
            for ch in el:
                if _strip(ch.tag) == "v":
                    val = ch.text
                elif _strip(ch.tag) == "is":
                    val = "".join(x.text or "" for x in ch.iter() if _strip(x.tag) == "t")
            if t == "s" and val is not None:
                val = ss[int(val)]
            elif t in (None, "n") and val is not None:
                try:
                    val = float(val)
                except ValueError:
                    pass
            r = int(re.search(r"\d+", ref).group(0)) - 1
            c = _col_idx(ref)
            rows.setdefault(r, {})[c] = val
            maxc = max(maxc, c + 1)
            el.clear()
        nrow = max(rows) + 1 if rows else 0
        out[name] = pd.DataFrame([[rows.get(i, {}).get(j) for j in range(maxc)] for i in range(nrow)])
    return out


# ---------------------------------------------------------------- helpers
rows, unresolved = [], []


def add_table(sig, claim_ids, kind, cell_type, cell_state, comparison, df, gene_col, stat_col,
              padj_col, source, selection):
    df = df.reset_index(drop=True)
    for i, r in df.iterrows():
        rows.append(dict(signature_id=f"{PAPER}__{sig}", paper=PAPER, claim_ids=";".join(claim_ids),
                         kind=kind, cell_type=cell_type, cell_state=cell_state, comparison=comparison,
                         gene=str(r[gene_col]).strip(), rank=i + 1,
                         stat=r[stat_col] if stat_col else "", padj=r[padj_col] if padj_col else "",
                         source=source, selection=selection))


def add_text(sig, claim_ids, kind, cell_type, cell_state, comparison, genes, selection="text list"):
    add_table(sig, claim_ids, kind, cell_type, cell_state, comparison,
              pd.DataFrame({"gene": genes}), "gene", None, None, "text", selection)


def unres(claim_id, term, reason):
    unresolved.append(dict(claim_id=claim_id, term=term, reason=reason))


# ---------------------------------------------------------------- DA subtype markers (Supp. Table 8)
book = read_strict_xlsx(XLSX, sheets={"Supplementary_Table_8"})
t8 = book["Supplementary_Table_8"]
t8.columns = [str(c) for c in t8.iloc[0]]
t8 = t8.iloc[1:].rename(columns={"None": "idx"})
t8["coef"] = t8["coef"].astype(float)
t8["fdr"] = t8["fdr"].astype(float)

extra_claims = {
    "SOX6_AGTR1": ["kamath2022-C06", "kamath2022-C09", "kamath2022-C12", "kamath2022-C17", "kamath2022-C18"],
    "CALB1_GEM": ["kamath2022-C07", "kamath2022-C13"],
    "CALB1_TRHR": ["kamath2022-C08", "kamath2022-C13"],
    "CALB1_RBP4": ["kamath2022-C13"],
}
SEL = "fdr<0.05 & coef>0 (MAST, subtype vs other DA subtypes), top 100 by coef"
for sub, d in t8.groupby("DA_subtype", sort=True):
    d = d[(d.fdr < 0.05) & (d.coef > 0)].sort_values("coef", ascending=False).head(100)
    add_table(f"DA_{sub}_markers", ["kamath2022-C01"] + extra_claims.get(sub, []), "state_markers",
              "DA neurons", sub, "", d, "primerid", "coef", "fdr",
              f"{XLSX.name}|Supplementary_Table_8", SEL)

# ---------------------------------------------------------------- text-only signatures
add_text("MG_GPNMB_textmarkers", ["kamath2022-C10"], "state_markers", "Microglia",
         "GPNMB+ (GPNMB/SULT1C2) microglia", "", ["GPNMB", "SULT1C2"])
unres("kamath2022-C10", "MG_GPNMB_SULT1C2 marker table",
      "marker table not in supplement (only DA subtype markers given); state-name genes used")
add_text("Astro_VIM_LHX2_textmarkers", ["kamath2022-C11"], "state_markers", "Astrocytes", "VIM_LHX2", "",
         ["VIM", "LHX2"])
unres("kamath2022-C11", "Astro_VIM_LHX2 marker table",
      "marker table not in supplement (only DA subtype markers given); state-name genes used")
add_text("CALB1_TMEM200A_smFISH_named", ["kamath2022-C13"], "named_genes", "DA neurons",
         "CALB1+/TMEM200A+", "", ["TH", "CALB1", "TMEM200A"])
add_text("SOX6_AGTR1_GWAS_named", ["kamath2022-C18"], "named_genes", "DA neurons", "SOX6_AGTR1", "",
         ["SNCA", "MAPT", "GAK", "WNT3", "IGSF9B"])
add_text("SOX6_AGTR1_TP53_named", ["kamath2022-C22"], "named_genes", "DA neurons", "SOX6_AGTR1",
         "PD/LBD vs control", ["TP53"], "text list (claimed up in PD/LBD)")
unres("kamath2022-C22", "SOX6_AGTR1 PD-vs-control DE table",
      "DE table not in supplement (Supp. Table 10 gives DEG counts only); TP53 from text used")

# TF-target gene sets: only TF names given
for cid, tf, d in [("kamath2022-C19", "TP53", "up"), ("kamath2022-C20", "NR2F2", "up"),
                   ("kamath2022-C21", "LMX1A", "down")]:
    unres(cid, f"{tf} targets (enrichR TRRUST/ENCODE/ARCHS4 TF target sets; {d} in SOX6_AGTR1)",
          "only TF name given; resolve target set centrally (TRRUST / ENCODE / ARCHS4 libraries)")


# ---------------------------------------------------------------- write + checks
def write_and_check():
    OUT.mkdir(parents=True, exist_ok=True)
    sig = pd.DataFrame(rows, columns=COLS)
    un = pd.DataFrame(unresolved, columns=["claim_id", "term", "reason"])
    assert sig.gene.str.len().gt(0).all() and not sig.gene.isin(["nan", "None"]).any(), "empty gene"
    assert (sig.groupby("signature_id")["rank"].min() == 1).all(), "rank must start at 1"
    sig.to_csv(OUT / f"{PAPER}.tsv", sep="\t", index=False)
    un.to_csv(OUT / f"{PAPER}_unresolved.tsv", sep="\t", index=False)

    print(sig.groupby(["signature_id", "kind"]).size().rename("n_genes").to_string())
    print("\nsignatures by kind:", sig.drop_duplicates("signature_id").kind.value_counts().to_dict())
    datelike = sig.gene[sig.gene.str.match(r"^\d{1,2}-[A-Za-z]{3}$|^[A-Za-z]{3}-\d{1,2}$|^\d{4}-\d{2}-\d{2}")]
    print("date-like gene symbols:", sorted(set(datelike)) or "none")

    claims = yaml.safe_load(open(ROOT / "claims" / f"{PAPER}.yaml"))
    need = [c["id"] for c in claims if c["test_template"] in ("state_signature", "geneset_shift", "de_pseudobulk")
            or (c["test_template"].startswith("composition") and (c.get("cell_state") or "").strip())]
    in_sig = set(";".join(sig.claim_ids).split(";"))
    in_un = set(un.claim_id)
    missing = [c for c in need if c not in in_sig and c not in in_un]
    print(f"\ncoverage: {len(need)} required claims; {len([c for c in need if c in in_sig])} in signatures; "
          f"{len([c for c in need if c not in in_sig and c in in_un])} unresolved only; missing: {missing or 'none'}")
    assert not missing


write_and_check()

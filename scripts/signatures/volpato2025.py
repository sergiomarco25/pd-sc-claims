"""Gene signatures for volpato2025 claims (see claims/SIGNATURES.md).

The medRxiv supplementary tables (Tables S2-S9: subtype markers, DEGs, pathway results) could not be
downloaded (medRxiv blocks scripts). The user supplied raw/supplements/volpato2025/
volpato2025_supplement_media-1.pdf, which holds supplementary FIGURES only (no gene tables); the few
gene names printed in those figures are used where they support a claim. Everything else comes from
the gene names in the claim text (named_genes, or state_markers with selection "text list"), and the
rest is listed as unresolved.

If xlsx/zip files are later dropped into the supplement folder, the script lists their sheets and
stops with a message, so the table parsing can be added (it cannot be written blind).

GEO metadata (GSE329625_full_metadata.csv.gz) is read only to confirm that the authors' subtype
labels (DaN_0-3, ODC_0-2, Astrocyte_0-3, Microglia_0-2) exist per cell, for later label-based re-tests.

Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/volpato2025.py
"""
import subprocess
import sys
import zipfile
from pathlib import Path

import pandas as pd
import yaml

PAPER = "volpato2025"
ROOT = Path(__file__).resolve().parents[2]

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

RAW = cfg.RAW
SUPP = RAW / "supplements" / PAPER
PDF = SUPP / "volpato2025_supplement_media-1.pdf"
GEO_META = RAW / PAPER / "GSE329625_full_metadata.csv.gz"
PDF_PY = Path.home() / ".local/share/mamba/envs/deck_rev/bin/python"   # has pymupdf
OUT = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison",
        "gene", "rank", "stat", "padj", "source", "selection"]
NO_SUPP = "supplement tables not available (medRxiv blocks download; supplied PDF has figures only)"

# ---------------------------------------------------------------- supplement tables present?
tables = sorted(p for p in SUPP.rglob("*") if p.suffix.lower() in (".xlsx", ".xls", ".zip", ".csv", ".tsv"))
if tables:
    print("Supplement table files found in", SUPP)
    for p in tables:
        if p.suffix.lower() == ".zip":
            print(" ", p.name, "->", zipfile.ZipFile(p).namelist())
        elif p.suffix.lower() in (".xlsx", ".xls"):
            print(" ", p.name, "->", pd.ExcelFile(p).sheet_names)
        else:
            print(" ", p.name)
    sys.exit("New volpato2025 supplement tables found: extend this script to parse them "
             "(extract zips into a subfolder of the paper folder), then rerun.")

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


def add_text(sig, claim_ids, kind, cell_type, cell_state, comparison, genes, selection="text list", source="text"):
    add_table(sig, claim_ids, kind, cell_type, cell_state, comparison,
              pd.DataFrame({"gene": genes}), "gene", None, None, source, selection)


def unres(claim_id, term, reason):
    unresolved.append(dict(claim_id=claim_id, term=term, reason=reason))


C = lambda *n: [f"{PAPER}-C{i:02d}" for i in n]

# ---------------------------------------------------------------- figure text from the supplied PDF
PDF_CODE = "import sys, fitz\nd = fitz.open(sys.argv[1])\nsys.stdout.write('\\n'.join(p.get_text() for p in d))"
pdf_text = ""
if PDF.exists():
    try:
        import fitz  # noqa: F401
        py = sys.executable
    except ImportError:
        py = str(PDF_PY)
    pdf_text = subprocess.run([py, "-I", "-c", PDF_CODE, str(PDF)], capture_output=True, text=True,
                              check=True).stdout
pdf_tokens = set(pdf_text.split())


def from_pdf(genes, fig):
    """Keep the figure-printed genes that are actually present in the PDF text."""
    found = [g for g in genes if g in pdf_tokens]
    missing = sorted(set(genes) - set(found))
    if missing:
        print(f"WARNING: genes not found in PDF text ({fig}):", missing)
    return found


# ---------------------------------------------------------------- GEO labels (check only)
if GEO_META.exists():
    meta = pd.read_csv(GEO_META, usecols=["CellType", "CellSubType"])
    print("GEO CellSubType labels:", sorted(meta.CellSubType.unique()))

# ---------------------------------------------------------------- subtype markers named in text (C01 etc.)
SUB = [("DaN_0_textmarkers", C(1), "DA neurons", "DaN_0", ["GRIA3", "DCX"]),
       ("DaN_1_textmarkers", C(1), "DA neurons", "DaN_1", ["SLC6A3", "GFRA2", "RET", "PITX3"]),
       ("DaN_2_textmarkers", C(1, 12), "DA neurons", "DaN_2", ["TMEFF2"]),
       ("DaN_3_textmarkers", C(1, 6, 9, 20, 24, 25, 26), "DA neurons", "AGTR1+ DaN_3",
        ["AGTR1", "SLC6A3", "GFRA2", "RET", "PITX3"]),
       ("psODC_2_textmarkers", C(1, 8), "Oligodendrocytes", "psODC_2 (ODC_2)",
        ["RBFOX1", "AFF3", "FMN1", "PDGFRA", "PDGFRB", "PDGFC", "CLDN10"]),
       ("Astrocyte_2_textmarkers", C(7), "Astrocytes", "Astrocyte_2", ["PTPRT", "ADGRV1", "EPHB1"]),
       ("Astrocyte_1_textmarkers", C(14, 26), "Astrocytes", "Astrocyte_1", ["GABRA2", "EDNRB", "SLC6A11"]),
       ("ODC_0_textmarkers", C(15), "Oligodendrocytes", "ODC_0", ["PLXDC2", "PALM2", "LAMA2"])]
for sig, cids, ct, state, genes in SUB:
    add_text(sig, cids, "state_markers", ct, state, "", genes,
             "text list (DaN_3 = DaN_1 markers + AGTR1-high, per text)" if state == "AGTR1+ DaN_3" else "text list")

add_text("C01_other_named", C(1), "named_genes", "All", "DaN / ODC / microglia subtypes", "",
         ["LMX1A", "LMX1B", "GRIK3", "OPALIN", "PLP1", "GPNMB", "SPP1", "TREM2"],
         "text list (genes named in Fig. 1e,f text not assigned to one subtype above)")
add_text("AGTR1_named", C(9), "named_genes", "DA neurons", "AGTR1+", "PD (Braak 5-6) vs control", ["AGTR1"],
         "text list (down in PD by RNAscope/IHC)")
add_text("DaN_3_dopamine_RAS_named", C(20), "named_genes", "DA neurons", "AGTR1+ DaN_3", "", ["AGTR1", "TH", "SLC6A3"])
add_text("DaN_3_T2D_targets_named", C(24), "named_genes", "DA neurons", "AGTR1+ DaN_3", "",
         ["GLP1R", "GCGR", "GIPR"], "text list (specific to DaN_3)")
add_text("PD_T2D_loci_named", C(25), "named_genes", "DA neurons", "AGTR1+ DaN_3", "", ["AGTR1", "TCF7L2"])
tcf = from_pdf(["TCF7L2"], "Fig. S16/S17")
if tcf:
    add_text("TCF7L2_regulon_named", C(26), "named_genes", "DA neurons", "AGTR1+ DaN_3; Astrocyte_1", "", tcf,
             "gene printed in Fig. S16 (expression) / S17 (SCENIC regulon AUC)", source=f"{PDF.name}|Fig. S16-S17")
else:
    add_text("TCF7L2_regulon_named", C(26), "named_genes", "DA neurons", "AGTR1+ DaN_3; Astrocyte_1", "",
             ["TCF7L2"])
opc = from_pdf(["PDGFRA", "PCDH15", "PTPRZ1", "VCAN", "BMP4", "SCN2A"], "Fig. S9")
if opc:
    add_text("OPC_differentiation_stage_named", C(17), "named_genes", "OPCs", "OPC -> psODC_2 / ODC_1", "", opc,
             "genes printed in Fig. S9 (OPC differentiation-stage markers)", source=f"{PDF.name}|Fig. S9")

# ---------------------------------------------------------------- unresolved
unres("volpato2025-C01", "full subtype marker lists (Table S2)", NO_SUPP + "; text-named markers used")
for cid, state in [("volpato2025-C06", "DaN_3"), ("volpato2025-C07", "Astrocyte_2"), ("volpato2025-C08", "psODC_2")]:
    unres(cid, f"{state} marker list (Table S2)", NO_SUPP + "; text-named markers used")
unres("volpato2025-C12", "per-subtype PD/iLBD-vs-control DEGs (Table S4)", NO_SUPP)
unres("volpato2025-C13", "DaN_1 PD DEGs (Table S4) and GO terms: cell-cell signalling, cell communication, "
      "mitochondrial processes / oxidative phosphorylation (Table S6)", NO_SUPP + "; resolve GO/Hallmark terms centrally")
unres("volpato2025-C14", "GO response to cytokine (GO:0034097) among Astrocyte_1 PD DEGs (Table S4)",
      NO_SUPP + "; resolve GO term centrally")
unres("volpato2025-C15", "GO response to unfolded protein (GO:0006986), cell activation, regulation of cell death "
      "among ODC_0 PD DEGs; stress response signalling in ODC_1 (Table S4)", NO_SUPP + "; resolve GO terms centrally")
unres("volpato2025-C24", "DaN_3 specificity statistics for GLP1R/GCGR/GIPR (Table S8)", NO_SUPP)
unres("volpato2025-C26", "SCENIC regulon target genes of TCF7L2 (Table S9)", NO_SUPP + "; regulon members not given")
unres("volpato2025-C27", "KEGG PI3K-Akt signalling pathway (hsa04151) / insulin resistance in perineuronal ODCs",
      "claim from medRxiv v2 abstract only; v2 full text and tables not accessed; resolve KEGG term centrally")


# ---------------------------------------------------------------- write + checks
def write_and_check():
    OUT.mkdir(parents=True, exist_ok=True)
    sig = pd.DataFrame(rows, columns=COLS)
    un = pd.DataFrame(unresolved, columns=["claim_id", "term", "reason"])
    assert sig.gene.str.len().gt(0).all() and not sig.gene.isin(["nan", "None"]).any(), "empty gene"
    assert (sig.groupby("signature_id")["rank"].min() == 1).all(), "rank must start at 1"
    assert not sig.duplicated(["signature_id", "gene"]).any(), "duplicate gene in a signature"
    sig.to_csv(OUT / f"{PAPER}.tsv", sep="\t", index=False)
    un.to_csv(OUT / f"{PAPER}_unresolved.tsv", sep="\t", index=False)

    print(sig.groupby(["signature_id", "kind"]).size().rename("n_genes").to_string())
    print("\nsignatures by kind:", sig.drop_duplicates("signature_id").kind.value_counts().to_dict())

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

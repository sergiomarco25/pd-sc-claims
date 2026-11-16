"""Gene signatures for smajic2022 claims (see claims/SIGNATURES.md).

Sources (raw/supplements/smajic2022/):
  awab446_supplementary_data.zip -> extracted/
      "awab446 Supplementary_data 1.pdf"  Supplementary Table 3 (cell-type markers, Seurat ROC test:
                                          myAUC, avg_diff, power, pct.1, pct.2; no p-values) as PDF text
      "awab446 Supplementary data 2.xlsx" Supplementary Table 5 (MAGMA per-gene ZSTAT)
  Supplementary_Table_8.xlsx.zip         Supplementary Table 8, which the PDF says is hosted on the
                                          authors' GitHub (SpielmannLab/pd_human_midbrain_snrnaseq,
                                          scrnaseq/data/Supplementary Table 8.xlsx.zip); downloaded
                                          2026-10-21 into the paper folder -> extracted/
      trajectory genes (Moran's I) with the "IPD differentially activated" flags (Fig. 3F, 4F, 4K, 4L),
      microglia GPNMBhigh / HSP90AA1high and astrocyte CD44high branch genes,
      per-cell-type IPD-vs-control DEGs (monocle3 GLM, q < 0.05) and GO enrichment of these DEGs.

The PDF is parsed with pypdf; if pypdf is not importable in this env, the text is extracted with
another local env that has it (PDF_PY below).

Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/smajic2022.py
"""
import re
import subprocess
import sys
import zipfile
from pathlib import Path

import pandas as pd
import yaml

PAPER = "smajic2022"
ROOT = Path(__file__).resolve().parents[2]

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

SUPP = cfg.RAW / "supplements" / PAPER
EXT = SUPP / "extracted"
PDF = EXT / "awab446 Supplementary_data 1.pdf"
T8 = EXT / "Supplementary Table 8.xlsx"
PDF_PY = Path.home() / "miniconda3/envs/troutpy_drvi/bin/python"
OUT = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison",
        "gene", "rank", "stat", "padj", "source", "selection"]

EXT.mkdir(exist_ok=True)
if not PDF.exists():
    zipfile.ZipFile(SUPP / "awab446_supplementary_data.zip").extractall(EXT)
if not T8.exists():
    zipfile.ZipFile(SUPP / "Supplementary_Table_8.xlsx.zip").extract(T8.name, EXT)

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


C = lambda *n: [f"{PAPER}-C{i:02d}" for i in n]

# ---------------------------------------------------------------- Supplementary Table 3 (PDF pages)
PDF_CODE = r"""
import sys, pypdf
r = pypdf.PdfReader(sys.argv[1])
for p in r.pages:
    sys.stdout.write((p.extract_text() or '') + '\n<<PAGE>>\n')
"""


def pdf_text(path):
    try:
        import pypdf  # noqa: F401
        return subprocess.run([sys.executable, "-I", "-c", PDF_CODE, str(path)], capture_output=True, text=True,
                              check=True).stdout
    except ImportError:
        return subprocess.run([str(PDF_PY), "-I", "-c", PDF_CODE, str(path)], capture_output=True, text=True,
                              check=True).stdout


def parse_table3(text):
    num = r"(-?[\d]+(?:[\.,]\d+)?(?:[eE]-?\d+)?)"
    full = re.compile(rf"^(.*?)\s+(\S+)\s+(ENSG\d+)\s+{num}\s+{num}\s+{num}\s+{num}\s+{num}\s*$")
    tail = re.compile(rf"^(ENSG\d+)\s+{num}\s+{num}\s+{num}\s+{num}\s+{num}\s*$")
    head = re.compile(r"^(.*?)\s+(\S+)\s*$")
    start = text.index("Supplementary Table 3 Gene expression markers")
    end = text.index("Supplementary Table 4 PD-associated", start)
    recs, buf = [], []
    for line in text[start:end].splitlines():
        line = line.strip()
        m = full.match(line)
        if m:
            recs.append(m.groups()); buf = []; continue
        m = tail.match(line)
        if m and buf:                          # gene row wrapped over lines (e.g. "ANKRD34C-" / "AS1")
            joined = ""
            for part in buf[-2:]:
                joined = joined + part if joined.endswith("-") else (joined + " " + part).strip()
            h = head.match(joined)
            recs.append((h.group(1), h.group(2)) + m.groups()); buf = []; continue
        if line and line != "<<PAGE>>":
            buf.append(line)
    df = pd.DataFrame(recs, columns=["cell_type", "gene", "ensembl", "myAUC", "avg_diff", "power", "pct1", "pct2"])
    for c in ["myAUC", "avg_diff", "power", "pct1", "pct2"]:
        df[c] = df[c].str.replace(",", ".").astype(float)
    return df


t3 = parse_table3(pdf_text(PDF))
print("Supp. Table 3 rows parsed per cell type:", t3.cell_type.value_counts().to_dict())
SRC3 = f"{PDF.name}|Supplementary Table 3"
SEL3 = "Seurat ROC markers listed by the authors (no p-values reported), avg_diff>0, top 100 by avg_diff"
for ct_pdf, ct, state, cids in [("CADPS2+ neurons", "Neurons (other)", "CADPS2high", C(1, 2)),
                                ("DaNs", "DA neurons", "", C(2)),
                                ("Microglia", "Microglia", "", C(23))]:
    d = t3[(t3.cell_type == ct_pdf) & (t3.avg_diff > 0)].sort_values("avg_diff", ascending=False).head(100)
    add_table(f"{ct_pdf.split()[0].replace('+', '')}_markers", cids, "state_markers", ct, state, "", d, "gene",
              "avg_diff", None, SRC3, SEL3)

# ---------------------------------------------------------------- Supplementary Table 8
x8 = pd.read_excel(T8, sheet_name=None, header=None)
SRC8 = T8.name + "|{}"


def gene_str(v):
    # Excel-converted symbols (e.g. MARCH1 -> 1-Mar) arrive as datetimes: keep them, as the date the cell holds
    return v.strftime("%Y-%m-%d") if hasattr(v, "strftime") else v


def table_at(sheet, header_row, cols, names):
    """Columns `cols` start with the gene-name column; the column left of it holds the Ensembl ID,
    which is used only when the gene name cell is empty."""
    raw = x8[sheet].iloc[header_row + 1:]
    d = raw.iloc[:, cols].copy()
    d.columns = names
    ens = raw.iloc[:, cols[0] - 1]
    d[names[0]] = d[names[0]].where(d[names[0]].notna(), ens.where(ens.astype(str).str.startswith("ENSG")))
    d[names[0]] = d[names[0]].map(gene_str)
    return d.dropna(subset=[names[0]])


# trajectory genes and branch genes
mt = table_at("Microglia_trajectory_genes", 3, [1, 2, 3, 4, 5], ["gene", "morans_test", "morans_I", "q", "ipd"])
mt["morans_I"] = mt.morans_I.astype(float)
mt_up = mt[mt.ipd.astype(str).str.strip() == "*"].sort_values("morans_I", ascending=False)
add_table("microglia_trajectory_IPD_up", C(10, 11), "de_up", "Microglia", "activation trajectory",
          "PD vs control (along trajectory)", mt_up, "gene", "morans_I", "q",
          SRC8.format("Microglia_trajectory_genes"),
          "trajectory genes flagged 'IPD differentially activated' (Fig. 3F), ranked by Moran's I")
for col, state in [(7, "GPNMBhigh microglia"), (9, "HSP90AA1high microglia")]:
    g = x8["Microglia_trajectory_genes"].iloc[4:, col].dropna().astype(str).str.strip()
    add_table(f"microglia_{state.split()[0]}_branch", C(10, 12), "state_markers", "Microglia", state, "",
              pd.DataFrame({"gene": g}), "gene", None, None, SRC8.format("Microglia_trajectory_genes"),
              "branch-activated genes listed by authors (Fig. 3E query genes), order as given")

at = table_at("Astrocyte_trajectory_genes", 2, [1, 2, 3, 4, 5], ["gene", "morans_test", "morans_I", "q", "ipd"])
at["morans_I"] = at.morans_I.astype(float)
at_up = at[at.ipd.astype(str).str.strip() == "*"].sort_values("morans_I", ascending=False)
add_table("astrocyte_trajectory_IPD_up", C(13, 15), "de_up", "Astrocytes", "CD44high",
          "PD vs control (along trajectory)", at_up, "gene", "morans_I", "q",
          SRC8.format("Astrocyte_trajectory_genes"),
          "trajectory genes flagged 'IPD differentially activated' (Fig. 4F), ranked by Moran's I")
g = x8["Astrocyte_trajectory_genes"].iloc[3:, 6].dropna().astype(str).str.strip()
add_table("astrocyte_CD44high_branch", C(13, 14, 15), "state_markers", "Astrocytes", "CD44high", "",
          pd.DataFrame({"gene": g}), "gene", None, None, SRC8.format("Astrocyte_trajectory_genes"),
          "branch-activated genes listed by authors (Fig. 4E query genes), order as given")

ot = table_at("Oligodendrocyte_trajectory_gene", 3, [1, 2, 3, 4, 5, 6], ["gene", "morans_test", "morans_I", "q", "up", "down"])
ot["morans_I"] = ot.morans_I.astype(float)
for kind, col, cids, fig in [("de_up", "up", C(16, 17, 18), "4K"), ("de_down", "down", C(16, 19), "4L")]:
    d = ot[ot[col].astype(str).str.strip() == "*"].sort_values("morans_I", ascending=False)
    add_table(f"oligodendrocyte_trajectory_IPD_{col}", cids, kind, "Oligodendrocytes", "",
              "PD vs control (along trajectory)", d, "gene", "morans_I", "q",
              SRC8.format("Oligodendrocyte_trajectory_gene"),
              f"trajectory genes flagged 'IPD differentially activated {col.upper()}' (Fig. {fig}), ranked by Moran's I")

# per-cell-type DEGs + GO terms
DE_SHEETS = {"Microglia_IPD_diff_exp_genes": ("Microglia", "Microglia", C(11)),
             "Astrocyte_IPD_diff_exp_genes": ("Astrocytes", "Astrocytes", C(14, 15)),
             "Oligodendrocytes_IPD_diff_exp_g": ("Oligodendrocytes", "Oligodendrocytes", C(17, 18, 19)),
             "Inhibitory_IPD_diff_exp_genes": ("Inhibitory", "Inhibitory neurons", []),
             "Excitatory_IPD_diff_exp_genes": ("Excitatory", "Excitatory neurons", []),
             "GABA_IPD_diff_exp_genes": ("GABA", "Inhibitory neurons", []),
             "DaNs_IPD_diff_exp_genes": ("DaNs", "DA neurons", []),
             "CADPS2high_IPD_diff_exp_genes": ("CADPS2high", "Neurons (other)", []),
             "OPC_IPD_diff_exp_genes": ("OPC", "OPCs", C(20)),
             "Ependymal_IPD_diff_exp_genes": ("Ependymal", "All", []),
             "Pericytes_IPD_diff_exp_genes": ("Pericytes", "Pericytes", C(20)),
             "Endothelial_cell_IPD_diff_exp_g": ("Endothelial", "Endothelial", C(20))}
SEL_DE = "all DEGs listed in Table 8 (monocle3 quasi-Poisson GLM on cells, q<0.05), split by sign of estimate, ranked by estimate"
go_all = []
for sheet, (short, ct, cids) in DE_SHEETS.items():
    raw = x8[sheet]
    if sheet == "Oligodendrocytes_IPD_diff_exp_g":
        de = pd.concat([table_at(sheet, 3, [1, 2, 3], ["gene", "q", "estimate"]),
                        table_at(sheet, 3, [6, 7, 8], ["gene", "q", "estimate"])])
        go_cols = [12, 15, 14]
    else:
        de = table_at(sheet, 3, [1, 2, 3], ["gene", "q", "estimate"])
        go_cols = [7, 10, 9]
    de["q"] = de.q.astype(float); de["estimate"] = de.estimate.astype(float)
    state = "CADPS2high" if short == "CADPS2high" else ("Ependymal" if short == "Ependymal" else "")
    for kind, sign in (("de_up", 1), ("de_down", -1)):
        d = de[(de.q < 0.05) & (sign * de.estimate > 0)].sort_values("estimate", ascending=(sign < 0))
        if d.empty:
            continue
        add_table(f"{short}_IPDvsCtrl_{kind[3:]}", cids, kind, ct, state, "PD vs control", d, "gene", "estimate",
                  "q", SRC8.format(sheet), SEL_DE)
    go = raw.iloc[4:, go_cols].dropna()
    go.columns = ["desc", "set", "genes"]
    go["sheet"], go["ct"] = sheet, ct
    go_all.append(go)
go_all = pd.concat(go_all)


def add_go(sig, cids, cell_type, sets, regex, desc):
    d = go_all[go_all.set.isin(sets) & go_all.desc.str.contains(regex, case=False)]
    genes = d.genes.str.split(",").explode().str.strip()
    cnt = genes[genes != ""].value_counts().rename_axis("gene").reset_index(name="n")
    cnt = cnt.sort_values(["n", "gene"], ascending=[False, True])
    add_table(sig, cids, "pathway_genes", cell_type, "", "PD vs control", cnt, "gene", "n", None,
              SRC8.format(";".join(sorted(d.sheet.unique()))),
              f"{desc}: union of genes of GO terms [{'; '.join(sorted(d.desc.unique()))}] in "
              f"[{', '.join(sorted(d.set.unique()))}] (GO enrichment of the IPD DEGs, p.adj<0.1); stat = n terms")


UPR = r"unfolded protein"
add_go("astrocyte_UPR_up", C(15), "Astrocytes", ["Astrocytes_up"], UPR, "UPR terms, astrocyte up-DEGs")
add_go("oligodendrocyte_UPR_up", C(18), "Oligodendrocytes", ["Oligodendrocytes_up"], UPR,
       "UPR terms, oligodendrocyte up-DEGs (trajectory-gene GO in Fig. 4K not tabulated)")
add_go("OPC_UPR_up", C(20), "OPCs", ["OPCs_up"], UPR, "UPR terms, OPC up-DEGs")
add_go("vascular_UPR_up", C(20), "Pericytes", ["Pericytes_up", "Endothelial_cells_up"], UPR,
       "UPR terms, pericyte + endothelial up-DEGs")
add_go("microglia_UPR_up", C(12), "Microglia", ["Microglia_up"], UPR,
       "UPR terms, microglia up-DEGs (branch-gene GO of Fig. 3E not tabulated)")

# ---------------------------------------------------------------- named genes (text)
add_text("CADPS2high_named", C(1, 2), "named_genes", "Neurons (other)", "CADPS2high", "",
         ["CADPS2", "TH", "TIAM1", "MAP2", "SCN2A"], "text list (CADPS2 high, TH low, TIAM1 higher)")
add_text("microglia_activation_named", C(11), "named_genes", "Microglia", "", "PD vs control",
         ["IL1B", "GPNMB", "HSP90AA1"], "text list (up in PD)")
add_text("microglia_branch_named", C(12), "named_genes", "Microglia", "GPNMBhigh / HSP90AA1high microglia", "",
         ["GPNMB", "HSP90AA1"])
add_text("astrocyte_CD44_named", C(14), "named_genes", "Astrocytes", "CD44high", "PD vs control", ["CD44"],
         "text list (up in PD)")
add_text("oligodendrocyte_states_named", C(16), "named_genes", "Oligodendrocytes", "OPALINhigh", "",
         ["OPALIN", "FRY", "RBFOX1", "S100B"], "text list (OPALIN/FRY myelinating state down; RBFOX1/S100B stress state up)")
add_text("oligodendrocyte_S100B_named", C(17), "named_genes", "Oligodendrocytes", "S100Bhigh", "PD vs control",
         ["S100B"], "text list (up in PD)")
add_text("GWAS_driver_named", C(23), "named_genes", "Microglia", "", "", ["LRRK2", "SNCA"],
         "text list (LRRK2 top in microglia/OPCs, SNCA in DaNs)")

# ---------------------------------------------------------------- unresolved
unres("smajic2022-C16", "OPALINhigh / RBFOX1-S100Bhigh oligodendrocyte subcluster markers",
      "subcluster marker table not in supplement (Fig. 4H only); named genes + trajectory genes used")
unres("smajic2022-C19", "GO neuron projection / synaptic transmission (Fig. 4L, oligodendrocyte trajectory-down genes)",
      "GO result for trajectory genes not tabulated and no such term among GO of oligodendrocyte down-DEGs; "
      "resolve GO terms centrally (trajectory-down gene list given as de_down)")
unres("smajic2022-C25", "astrocyte / microglia proliferation",
      "no genes or supporting Results given (abstract only); a proliferation set (e.g. GO cell cycle) must be chosen centrally")
unres("smajic2022-C12", "GO enrichment of microglia branch genes (Fig. 3E: cytokine secretion, UPR)",
      "branch-gene GO table not in supplement; branch gene lists + UPR terms of microglia DEGs used")


# ---------------------------------------------------------------- write + checks
def write_and_check():
    OUT.mkdir(parents=True, exist_ok=True)
    sig = pd.DataFrame(rows, columns=COLS)
    un = pd.DataFrame(unresolved, columns=["claim_id", "term", "reason"])
    assert sig.gene.str.len().gt(0).all() and not sig.gene.isin(["nan", "None"]).any(), "empty gene"
    assert (sig.groupby("signature_id")["rank"].min() == 1).all(), "rank must start at 1"
    dup = sig[sig.duplicated(["signature_id", "gene"], keep=False)]
    print("duplicate genes within a signature (kept):", dup[["signature_id", "gene"]].values.tolist() or "none")
    sig.to_csv(OUT / f"{PAPER}.tsv", sep="\t", index=False)
    un.to_csv(OUT / f"{PAPER}_unresolved.tsv", sep="\t", index=False)

    print(sig.groupby(["signature_id", "kind"]).size().rename("n_genes").to_string())
    print("\nsignatures by kind:", sig.drop_duplicates("signature_id").kind.value_counts().to_dict())
    datelike = sig[sig.gene.str.match(r"^\d{1,2}-[A-Za-z]{3}$|^[A-Za-z]{3}-\d{1,2}$|^\d{4}-\d{2}-\d{2}|^\d{5}$")]
    print("date-like gene symbols (kept as reported):", datelike[["signature_id", "gene"]].values.tolist() or "none")

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

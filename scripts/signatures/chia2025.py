"""Signatures for chia2025 (Chia et al. 2025, Genome Med; PFC snRNA, 4 NCI / 4 AD / 6 DLB / 6 PDD).

Reads 13073_2025_1519_MOESM1_ESM.xlsx (Additional file 1, Tables S1-S10) and writes
claims/signatures/chia2025.tsv and claims/signatures/chia2025_unresolved.tsv.
Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/chia2025.py

Notes
- DEG tables S5-S9 list the union of DESeq2 pseudobulk (nominal P<0.05) and MAST (BH) hits with a
  `keep` flag. keep == True is the authors' final DEG set: it reproduces the reported 294
  DLB-vs-PDD DEGs with 162 in microglia (claim C13). Only keep == True rows are used.
- stat = avg_log2FC (authors' combined column: mean of bulk and MAST log2FC where both exist);
  padj = p_adjust_MAST when present (empty for pseudobulk-only hits, whose BH-adjusted p is ~1).
- Table S9's title row is a copy of S8's ("LBD vs NCI"); the sheet name says DLB vs PDD and the
  counts match the DLB-vs-PDD text, so it is treated as DLB vs PDD. The reference group of that
  comparison is not stated; sign assumed positive = higher in DLB.
- No subcluster marker tables (SST/SV2C interneurons, 'distinct microglia') and no pathway tables
  are in the supplement; those claims get named genes from the text + unresolved rows.
"""
from pathlib import Path
import pandas as pd
import yaml

PAPER = "chia2025"
ROOT = Path(__file__).resolve().parents[2]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

SUP = cfg.RAW / "supplements" / PAPER
XL = SUP / "13073_2025_1519_MOESM1_ESM.xlsx"
OUT = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison",
        "gene", "rank", "stat", "padj", "source", "selection"]
C = lambda *n: [f"{PAPER}-C{int(i):02d}" for i in n]  # noqa: E731
CT = {"Astrocyte": ("Astrocytes", "Astro"), "Excitatory neuron": ("Excitatory neurons", "ExN"),
      "Inhibitory neuron": ("Inhibitory neurons", "InN"), "Microglia": ("Microglia", "Micro"),
      "OPC": ("OPCs", "OPC"), "Oligodendrocyte": ("Oligodendrocytes", "Oligo"),
      "Endothelial cell": ("Endothelial", "Endo")}

rows, unresolved, DATE_GENES = [], [], set()


def gname(g):
    if hasattr(g, "strftime"):  # Excel date-corrupted symbol: keep as Excel displays it, flag
        s = f"{g.day}-{g.strftime('%b')}"
        DATE_GENES.add(s)
        return s
    return str(g).strip()


def add(name, claims, kind, ct, state, comp, genes, stats, padjs, source, selection):
    genes = [gname(g) for g in genes]
    for i, g in enumerate(genes):
        rows.append(dict(signature_id=f"{PAPER}__{name}", paper=PAPER, claim_ids=";".join(claims),
                         kind=kind, cell_type=ct, cell_state=state, comparison=comp, gene=g, rank=i + 1,
                         stat="" if stats is None else stats[i], padj="" if padjs is None else padjs[i],
                         source=source, selection=selection))


def unres(claim, term, reason):
    unresolved.append(dict(claim_id=claim, term=term, reason=reason))


def table(sheet):
    raw = pd.read_excel(XL, sheet_name=sheet, header=None)
    h = raw.index[raw.iloc[:, 1].astype(str).str.strip() == "gene"][0]
    d = raw.iloc[h + 1:, 1:].copy()
    d.columns = [str(c).strip() for c in raw.iloc[h, 1:]]
    d = d[d.gene.notna() & d.cell_type.notna()].reset_index(drop=True)
    for c in ("avg_log2FC", "p_adjust_MAST"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    return d


# sheet -> (short comparison, label, base claims, extra claims per cell type)
DE = {
    "Table S7 PDD DEG": ("PDD_vs_NCI", "PDD vs control", C(5, 6),
                         {"Oligodendrocyte": C(7, 8), "Microglia": C(8, 11, 12)}),
    "Table S6 DLB DEG": ("DLB_vs_NCI", "DLB vs control", C(5),
                         {"Oligodendrocyte": C(7, 8), "Microglia": C(8, 11, 12)}),
    "Table S5 AD DEG": ("AD_vs_NCI", "AD vs control", C(5),
                        {"Oligodendrocyte": C(7, 8), "Microglia": C(8, 11)}),
    "Table S8 DLB PDD vs NCI": ("LBD_vs_NCI", "LBD (PDD + DLB) vs control", C(5), {}),
    "Table S9 DLB vs PDD": ("DLB_vs_PDD", "DLB vs PDD", C(5, 13), {}),
}
for sheet, (short, lbl, base, extra) in DE.items():
    d = table(sheet)
    d = d[d.keep.astype(str).str.strip().str.lower() == "true"]
    print(f"{sheet}: {len(d)} keep==True DEGs; per cell type {d.cell_type.value_counts().to_dict()}")
    for ct, dc in d.groupby("cell_type"):
        ctn, cts = CT[ct]
        cl = sorted(set(base + extra.get(ct, [])))
        for kind, sub in (("de_up", dc[dc.avg_log2FC > 0]), ("de_down", dc[dc.avg_log2FC < 0])):
            sub = sub.reindex(sub.avg_log2FC.abs().sort_values(ascending=False).index)
            if len(sub):
                add(f"{cts}_{short}_{kind.split('_')[1]}", cl, kind, ctn, "", lbl, sub.gene,
                    sub.avg_log2FC.round(4).tolist(),
                    ["" if pd.isna(p) else p for p in sub.p_adjust_MAST], f"{XL.name}|{sheet.strip()}",
                    "authors' final DEG set (keep==True: DESeq2 pseudobulk P<0.05 union MAST BH); "
                    "ranked by |avg_log2FC|; padj = MAST BH (empty if pseudobulk-only)")

# ---------------- named genes ----------------
named = [
    ("InN_SST_SV2C_named", C(4), "Inhibitory neurons", "SST IN, SV2C IN", "PDD vs control; DLB vs control", ["SST", "SV2C"]),
    ("Oligo_HSPA1A_up", C(7), "Oligodendrocytes", "", "PDD vs control", ["HSPA1A"]),
    ("distinct_microglia_named", C(9), "Microglia", "distinct microglia (microglia/oligodendrocyte marker co-expressing)", "",
     ["P2RY12", "PTPRC", "APOE", "CD74", "PLP1", "MBP", "ST18"]),
    ("Micro_MSR1_up", C(11), "Microglia", "", "PDD vs control", ["MSR1"]),
]
for name, cl, ct, st, comp, genes in named:
    add(name, cl, "named_genes", ct, st, comp, genes, None, None, "text", "text list")

# ---------------- unresolved ----------------
unres(C(4)[0], "SST IN and SV2C IN subcluster markers", "subcluster marker table not in supplement; only the "
      "defining genes SST and SV2C emitted")
unres(C(9)[0], "distinct microglia state markers", "marker table not in supplement; only text-named markers emitted "
      "(state may reflect microglia-oligodendrocyte doublets / ambient myelin RNA)")
for t in ("cell differentiation / cell migration", "synapse development / synapse organisation",
          "lipid metabolism and lipid transport", "synaptic plasticity / learning or memory / cognition",
          "antigen processing and presentation"):
    unres(C(12)[0], f"GO BP: {t}", "pathway named only (Fig. 2e / 4g; no enrichment table in supplement); resolve "
          "from MSigDB; microglia DEG-up lists attached")
for t in ("synapse organisation", "neuron projection development", "cell adhesion", "ion transport"):
    unres(C(13)[0], f"GO BP: {t}", "pathway named only (Fig. 3b); resolve from MSigDB; DLB-vs-PDD DEG lists attached")

# ---------------- write + checks ----------------
OUT.mkdir(parents=True, exist_ok=True)
sig = pd.DataFrame(rows, columns=COLS)
un = pd.DataFrame(unresolved, columns=["claim_id", "term", "reason"])
assert sig.gene.str.len().gt(0).all() and not sig.gene.isin(["nan", "None"]).any()
assert (sig.groupby("signature_id")["rank"].min() == 1).all()
dup = sig.duplicated(["signature_id", "gene"])
if dup.any():
    print("WARNING duplicated genes within signature (kept):", sig[dup][["signature_id", "gene"]].values.tolist())
sig.to_csv(OUT / f"{PAPER}.tsv", sep="\t", index=False)
un.to_csv(OUT / f"{PAPER}_unresolved.tsv", sep="\t", index=False)

print(sig.groupby(["kind", "signature_id"]).size().to_string())
print("\nsignatures by kind:", sig.drop_duplicates("signature_id").kind.value_counts().to_dict())
susp = sig[sig.gene.str.match(r"^\d+-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$|^(Jan|Feb|Mar|Apr|Sep|Oct|Nov|Dec)-\d+$|^[\d.]+$")]
print("possible Excel date-corrupted / numeric symbols:", sorted(set(susp.gene) | DATE_GENES),
      sig[sig.gene.isin(set(susp.gene) | DATE_GENES)].signature_id.unique().tolist())

claims = yaml.safe_load(open(ROOT / "claims" / f"{PAPER}.yaml"))
need = {c["id"] for c in claims if c["test_template"] in ("state_signature", "geneset_shift", "de_pseudobulk")
        or (c["test_template"].startswith("composition") and (c.get("cell_state") or "").strip())}
in_sig = {x for s in sig.claim_ids.unique() for x in s.split(";")}
in_un = set(un.claim_id)
print(f"\ncoverage: {len(need)} claims need signatures; {len(need & in_sig)} in signatures; "
      f"{len(need & in_un - in_sig)} unresolved only; missing: {sorted(need - in_sig - in_un)}")

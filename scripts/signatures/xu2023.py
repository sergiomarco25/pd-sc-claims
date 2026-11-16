"""Signatures for xu2023 (Xu et al. 2023, Nat Aging; putamen snRNA, 4 CTRL / 4 AD / 4 PD).

Reads 43587_2023_363_MOESM3_ESM.xlsx (Supplementary Tables S1-S69) and writes
claims/signatures/xu2023.tsv and claims/signatures/xu2023_unresolved.tsv.
Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/xu2023.py

Notes
- Conserved-marker tables (FindConservedMarkers) report per-group p_val_adj and avg_log2FC for
  CTRL, PD and AD. A gene is kept when the largest of the three p_val_adj < 0.05 and all three
  log2FC > 0; ranked by mean log2FC across groups; padj column = max group p_val_adj.
- DEG tables (FindMarkers, Wilcoxon on cells) are already filtered by the authors at
  FDR < 0.05 & |avg_logFC| > 0.25 (natural log); positive avg_logFC = higher in disease
  (checked: CRYAB/HSPB1 up, APOE down in PD astrocytes, as in the text).
"""
from pathlib import Path
import re
import pandas as pd
import yaml

PAPER = "xu2023"
ROOT = Path(__file__).resolve().parents[2]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

SUP = cfg.RAW / "supplements" / PAPER
XL = SUP / "43587_2023_363_MOESM3_ESM.xlsx"
OUT = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison",
        "gene", "rank", "stat", "padj", "source", "selection"]
C = lambda *n: [f"{PAPER}-C{int(i):02d}" for i in n]  # noqa: E731

rows, unresolved = [], []
xl = pd.ExcelFile(XL)


def gname(g):
    """Gene symbol as displayed by Excel; date cells (Excel-corrupted symbols, e.g. MARCH1) become
    'd-Mon' as Excel shows them and are reported in DATE_GENES (kept, not renamed)."""
    if hasattr(g, "strftime"):
        s = f"{g.day}-{g.strftime('%b')}"
        DATE_GENES.add(s)
        return s
    return str(g).strip()


DATE_GENES = set()


def add(name, claims, kind, ct, state, comp, genes, stats, padjs, source, selection):
    genes = [gname(g) for g in genes]
    for i, g in enumerate(genes):
        rows.append(dict(signature_id=f"{PAPER}__{name}", paper=PAPER, claim_ids=";".join(claims),
                         kind=kind, cell_type=ct, cell_state=state, comparison=comp, gene=g, rank=i + 1,
                         stat="" if stats is None else stats[i], padj="" if padjs is None else padjs[i],
                         source=source, selection=selection))


def unres(claim, term, reason):
    unresolved.append(dict(claim_id=claim, term=term, reason=reason))


def sheet(prefix):
    """Return (sheet name, table with header row 'Gene')."""
    s = [s for s in xl.sheet_names if s.startswith(prefix + "_")]
    assert len(s) == 1, (prefix, s)
    raw = xl.parse(s[0], header=None)
    h = raw.index[raw.iloc[:, 0].astype(str).str.strip() == "Gene"][0]
    d = raw.iloc[h + 1:].copy()
    d.columns = [str(c).strip() for c in raw.iloc[h]]
    d = d[d.Gene.notna()].reset_index(drop=True)
    return s[0], d


def markers(tab, name, claims, ct, state):
    s, d = sheet(tab)
    groups = ["CTRL", "PD", "AD"]
    # S24 (Micr-1) labels the fold change avg_logFC (natural log); others avg_log2FC
    d = d.rename(columns={f"{g}_avg_logFC": f"{g}_avg_log2FC" for g in groups})
    for g in groups:
        for k in ("avg_log2FC", "p_val_adj"):
            d[f"{g}_{k}"] = pd.to_numeric(d[f"{g}_{k}"])
    d["padj"] = d[[f"{g}_p_val_adj" for g in groups]].max(axis=1)
    d["lfc_min"] = d[[f"{g}_avg_log2FC" for g in groups]].min(axis=1)
    d["lfc_mean"] = d[[f"{g}_avg_log2FC" for g in groups]].mean(axis=1)
    k = d[(d.padj < 0.05) & (d.lfc_min > 0)].sort_values("lfc_mean", ascending=False).head(100)
    add(name, claims, "state_markers", ct, state, "", k.Gene, k.lfc_mean.round(4).tolist(),
        k.padj.tolist(), f"{XL.name}|{s}",
        "conserved markers: max group p_val_adj<0.05 & log2FC>0 in CTRL, PD and AD; "
        "top 100 by mean log2FC (stat = mean log2FC, padj = max group p_val_adj)"
        + ("; S24 reports avg_logFC (natural log)" if tab == "Table_S24" else ""))
    return len(k)


def de(tab, name, claims, ct, state, comp):
    s, d = sheet(tab)
    d["avg_logFC"] = pd.to_numeric(d.avg_logFC)
    d["p_val_adj"] = pd.to_numeric(d.p_val_adj)
    d = d[(d.p_val_adj < 0.05) & (d.avg_logFC.abs() > 0.25)]
    for kind, sub in (("de_up", d[d.avg_logFC > 0]), ("de_down", d[d.avg_logFC < 0])):
        sub = sub.reindex(sub.avg_logFC.abs().sort_values(ascending=False).index)
        if len(sub):
            add(f"{name}_{kind.split('_')[1]}", claims, kind, ct, state, comp, sub.Gene,
                sub.avg_logFC.round(4).tolist(), sub.p_val_adj.tolist(), f"{XL.name}|{s}",
                "authors' DEG table: p_val_adj<0.05 & |avg_logFC|>0.25 (natural log, Wilcoxon on cells); "
                "ranked by |avg_logFC|")


# ---------------- astrocytes ----------------
for t, st, cl in (("Table_S2", "Ast-0", C(1)), ("Table_S3", "Ast-1", C(1, 2)), ("Table_S4", "Ast-2", C(1, 2))):
    markers(t, f"{st.replace('-', '')}_markers", cl, "Astrocytes", st)
for t, st, comp, cl in (("Table_S12", "Ast-0", "AD vs control", C(6)), ("Table_S13", "Ast-1", "AD vs control", C(6)),
                        ("Table_S14", "Ast-2", "AD vs control", C(6)),
                        ("Table_S15", "Ast-0", "PD vs control", C(3, 4, 6, 7, 8, 9)),
                        ("Table_S16", "Ast-1", "PD vs control", C(3, 4, 6, 7, 8, 9)),
                        ("Table_S17", "Ast-2", "PD vs control", C(3, 4, 6, 7, 8, 9))):
    de(t, f"{st.replace('-', '')}_{comp.split()[0]}_vs_CTRL", cl, "Astrocytes", st, comp)

# ---------------- microglia ----------------
markers("Table_S24", "Micr1_markers", C(11), "Microglia", "Micr-1")
for t, st, comp, cl in (("Table_S33", "Micr-0", "AD vs control", C(10)), ("Table_S34", "Micr-1", "AD vs control", C(10)),
                        ("Table_S35", "Micr-0", "PD vs control", C(10, 13, 14, 15)),
                        ("Table_S36", "Micr-1", "PD vs control", C(10, 13, 14, 15))):
    de(t, f"{st.replace('-', '')}_{comp.split()[0]}_vs_CTRL", cl, "Microglia", st, comp)

# ---------------- neurons ----------------
NEU = {"mD1": 38, "pD1": 39, "ncD1": 40, "mD2": 41, "pD2": 42, "ncD2": 43, "hMSN": 44, "sMSN": 45,
       "in.PDGFD": 46, "in.CALB2": 47, "in.TH": 48, "in.SST": 49, "in.SLC5A7": 50}
for st, n in NEU.items():
    ct = "Inhibitory neurons" if st.startswith("in.") else "MSNs"
    markers(f"Table_S{n}", f"{st.replace('.', '_')}_markers", C(17, 18) if st == "sMSN" else C(17), ct, st)
# DEG tables present in the workbook (S37 absent; no PD DEG table for ncD2, in.SST, in.SLC5A7)
present = {}
for s in xl.sheet_names:
    m = re.match(r"(Table_S(\d+))_DEG_(AD|PD)_vs_CTRL_(.+)$", s)
    if m and int(m.group(2)) >= 51:
        present[(m.group(3), m.group(4))] = m.group(1)
abbr = {"in.PDG": "in.PDGFD", "in.CAL": "in.CALB2"}
for (dis, st0), t in sorted(present.items()):
    st = abbr.get(st0, st0)
    ct = "Inhibitory neurons" if st.startswith("in.") else "MSNs"
    cl = C(19) + (C(20) if dis == "PD" and st in ("mD1", "pD1") else [])
    de(t, f"{st.replace('.', '_')}_{dis}_vs_CTRL", cl, ct, st, f"{dis} vs control")
missing_pd = [s for s in NEU if ("PD", s) not in {(d, abbr.get(x, x)) for d, x in present}]
unres(C(19)[0], "neuronal DEG tables PD vs control for " + ", ".join(missing_pd),
      "not in supplement (Supplementary Table S37 absent; other subclusters covered)")

# ---------------- named genes (claim text) ----------------
named = [
    ("Ast0_named", C(1), "Astrocytes", "Ast-0", "", ["GPC5", "NRXN1"]),
    ("Ast1_named", C(1), "Astrocytes", "Ast-1", "", ["APOE", "CLU", "APOC1"]),
    ("Ast2_named", C(1), "Astrocytes", "Ast-2", "", ["GFAP", "CD44", "TNC"]),
    ("Ast1_reactive_named", C(2), "Astrocytes", "Ast-1", "", ["S100B", "VIM", "MT2A", "MT1E", "CRYAB", "MT1G"]),
    ("Ast2_activated_named", C(2), "Astrocytes", "Ast-2", "", ["GFAP", "CD44", "C3", "SYNM", "MAOB"]),
    ("Astro_PD_up_named", C(3), "Astrocytes", "", "PD vs control", ["CRYAB", "HSPB1", "EGFR", "TLR4"]),
    ("Astro_PD_down_APOE", C(4), "Astrocytes", "", "PD vs control", ["APOE"]),
    ("Astro_PD_down_F3", C(5), "Astrocytes", "", "PD vs control", ["F3"]),
    ("Micr1_named", C(11), "Microglia", "Micr-1", "", ["AIF1", "APOC1", "CD14", "FTL", "APOE", "B2M", "TYROBP", "TREM2"]),
    ("Micr1_PD_up_TREM2", C(12), "Microglia", "Micr-1", "PD vs control", ["TREM2"]),
    ("Micro_PD_down_AB_binding", C(14), "Microglia", "", "PD vs control", ["MSR1", "CST3"]),
    ("Micro_PD_up_protein_folding", C(15), "Microglia", "", "PD vs control", ["HSP90AA1", "HSP90AB1", "FKBP4"]),
    ("MSN_subtype_named", C(17), "MSNs", "mD1/pD1/ncD1/mD2/pD2/ncD2/hMSN/sMSN", "",
     ["DRD1", "TAC1", "DRD2", "PENK", "PDYN", "RXFP2", "FAM210B", "MPP3", "POU6F2", "LGR5", "FHAD1",
      "RXFP1", "MKX", "JAG1", "STXBP6", "EPHA4", "KCNIP1", "BACH2"]),
    ("sMSN_named", C(18), "MSNs", "sMSN", "", ["HPCA", "PSAP", "DKK3"]),
]
for name, cl, ct, st, comp, genes in named:
    add(name, cl, "named_genes", ct, st, comp, genes, None, None, "text", "text list")

# ---------------- pathway terms without member genes ----------------
unres(C(7)[0], "KEGG neurodegeneration pathways (Parkinson / Alzheimer / Huntington disease; Fig. 3d)",
      "pathway named only (enrichment table not in supplement); resolve from MSigDB; astrocyte PD DEG-up lists attached")
unres(C(8)[0], "GO:0001540 amyloid-beta binding", "GO term named only; resolve from MSigDB; astrocyte PD DEG-down lists attached")
unres(C(9)[0], "AD/PD GWAS risk genes DE in >=2 astrocyte subclusters (Fig. 3e)",
      "gene list shown only in figure heatmap; PD GWAS gene list needed centrally; astrocyte DEG lists attached")
unres(C(14)[0], "GO:0001540 amyloid-beta binding", "GO term named only; resolve from MSigDB (MSR1, CST3 emitted as named genes)")
for term in ("GO:0060291 long-term synaptic potentiation", "GO:0048167 regulation of synaptic plasticity",
             "GO:0006836 neurotransmitter transport", "GO:0007611 learning or memory"):
    unres(C(20)[0], term, "GO term named only (Fig. 8i); resolve from MSigDB; mD1/pD1 PD DEG-down lists attached")

# ---------------- write + checks ----------------
OUT.mkdir(parents=True, exist_ok=True)
sig = pd.DataFrame(rows, columns=COLS)
un = pd.DataFrame(unresolved, columns=["claim_id", "term", "reason"])
assert sig.gene.str.len().gt(0).all() and not sig.gene.isin(["nan", "None"]).any()
assert (sig.groupby("signature_id")["rank"].min() == 1).all()
assert not sig.duplicated(["signature_id", "gene"]).any(), sig[sig.duplicated(["signature_id", "gene"], keep=False)]
sig.to_csv(OUT / f"{PAPER}.tsv", sep="\t", index=False)
un.to_csv(OUT / f"{PAPER}_unresolved.tsv", sep="\t", index=False)

print(sig.groupby(["kind", "signature_id"]).size().to_string())
print("\nsignatures by kind:", sig.drop_duplicates("signature_id").kind.value_counts().to_dict())
susp = sig[sig.gene.str.match(r"^\d+-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$|^(Jan|Feb|Mar|Apr|Sep|Oct|Nov|Dec)-\d+$|^\d{4}-\d{2}-\d{2}")]
print("possible Excel date-corrupted symbols:", sorted(set(susp.gene) | DATE_GENES),
      sig[sig.gene.isin(DATE_GENES)].signature_id.unique().tolist())

claims = yaml.safe_load(open(ROOT / "claims" / f"{PAPER}.yaml"))
need = {c["id"] for c in claims if c["test_template"] in ("state_signature", "geneset_shift", "de_pseudobulk")
        or (c["test_template"].startswith("composition") and (c.get("cell_state") or "").strip())}
in_sig = {x for s in sig.claim_ids.unique() for x in s.split(";")}
in_un = set(un.claim_id)
print(f"\ncoverage: {len(need)} claims need signatures; {len(need & in_sig)} in signatures; "
      f"{len(need & in_un - in_sig)} unresolved only; missing: {sorted(need - in_sig - in_un)}")

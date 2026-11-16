"""Signatures for dehestani2024 (Dehestani et al. 2024, Mol Brain; PFC + ACC snRNA,
2 GBA1-PD / 2 LRRK2-PD / 2 controls).

Reads
  13041_2024_1128_MOESM5_ESM.xlsx  Table S5: DEGs (Seurat FindMarkers, cells) ODC/OPC x GBA1/LRRK2 vs HC
  13041_2024_1128_MOESM7_ESM.xlsx  Table S7: Enrichr GO BP enrichment of up/down DEGs (Genes column)
and writes claims/signatures/dehestani2024.tsv and claims/signatures/dehestani2024_unresolved.tsv.
Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/dehestani2024.py

Notes
- Table S5 sheets are already filtered (p_val_adj < 0.05, |avg_log2FC| > 0.25). Positive avg_log2FC =
  higher in the mutation carriers (checked: the 'negative regulation of inclusion body assembly'
  genes of GBA_ODC_UP all have positive log2FC in ODC_GBA_HC).
- Pathway signatures take the Enrichr overlap genes of the GO term the claim names (matched by GO id).
  Several of these terms are not significant (Adjusted P-value >= 0.05) in Table S7; the value is
  written in `selection`.
"""
from pathlib import Path
import re
import pandas as pd
import yaml

PAPER = "dehestani2024"
ROOT = Path(__file__).resolve().parents[2]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

SUP = cfg.RAW / "supplements" / PAPER
X_DE = SUP / "13041_2024_1128_MOESM5_ESM.xlsx"
X_GO = SUP / "13041_2024_1128_MOESM7_ESM.xlsx"
OUT = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison",
        "gene", "rank", "stat", "padj", "source", "selection"]
C = lambda *n: [f"{PAPER}-C{int(i):02d}" for i in n]  # noqa: E731
CT = {"ODC": "Oligodendrocytes", "OPC": "OPCs"}
MUT = {"GBA": "GBA1", "LRRK2": "LRRK2"}
PW_CLAIM = {("ODC", "GBA"): C(5), ("OPC", "GBA"): C(6), ("ODC", "LRRK2"): C(7), ("OPC", "LRRK2"): C(8)}

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


# ---------------- Table S5: DEGs ----------------
lfc = {}
for ctk in ("ODC", "OPC"):
    for mut in ("GBA", "LRRK2"):
        sheet = f"{ctk}_{mut}_HC"
        d = pd.read_excel(X_DE, sheet_name=sheet)
        d = d[d.gene.notna() & (d.p_val_adj < 0.05) & (d.avg_log2FC.abs() > 0.25)]
        lfc[(ctk, mut)] = dict(zip(d.gene.map(gname), d.avg_log2FC))
        comp = f"{MUT[mut]}-PD vs control"
        cl = sorted(set(C(3) + PW_CLAIM[(ctk, mut)]))
        for kind, sub in (("de_up", d[d.avg_log2FC > 0]), ("de_down", d[d.avg_log2FC < 0])):
            sub = sub.reindex(sub.avg_log2FC.abs().sort_values(ascending=False).index)
            add(f"{ctk}_{MUT[mut]}_vs_HC_{kind.split('_')[1]}", cl, kind, CT[ctk], "", comp, sub.gene,
                sub.avg_log2FC.round(4).tolist(), sub.p_val_adj.tolist(), f"{X_DE.name}|{sheet}",
                "Table S5 (Wilcoxon on cells, PFC+ACC pooled): p_val_adj<0.05 & |avg_log2FC|>0.25; ranked by |avg_log2FC|")

# ---------------- Table S7: named GO terms ----------------
TERMS = [  # (sheet, GO id, claim)
    ("GBA_ODC_UP", "GO:0090084"), ("GBA_ODC_DOWN", "GO:0006627"),
    ("GBA_OPC_UP", "GO:1902949"), ("GBA_OPC_DOWN", "GO:0050804"),
    ("LRRK2_ODC_UP", "GO:0010968"), ("LRRK2_ODC_DOWN", "GO:0043266"),
    ("LRRK2_OPC_UP", "GO:0032271"), ("LRRK2_OPC_DOWN", "GO:0006695"),
]
for sheet, goid in TERMS:
    mut, ctk, dirn = sheet.split("_")
    t = pd.read_excel(X_GO, sheet_name=sheet)
    r = t[t.Term.astype(str).str.contains(re.escape(f"({goid})"))]
    assert len(r) == 1, (sheet, goid)
    r = r.iloc[0]
    desc = r.Term.split(" (GO:")[0]
    genes = [g.strip() for g in str(r.Genes).split(";") if g.strip()]
    lf = lfc[(ctk, mut)]
    genes = sorted(genes, key=lambda g: -abs(lf.get(g, 0)))
    claim = PW_CLAIM[(ctk, mut)]
    short = re.sub(r"[^A-Za-z0-9]+", "_", desc)[:45].strip("_")
    add(f"{ctk}_{MUT[mut]}_{dirn.lower()}_{short}", claim, "pathway_genes", CT[ctk], "",
        f"{MUT[mut]}-PD vs control ({dirn.lower()}regulated DEGs)", genes,
        [round(lf[g], 4) if g in lf else "" for g in genes], None, f"{X_GO.name}|{sheet}",
        f"Enrichr overlap genes of {goid} {desc} (overlap {r.Overlap}, adj. P {r['Adjusted P-value']:.3g}"
        f"{'' if r['Adjusted P-value'] < 0.05 else ', NOT significant'}); ranked by |avg_log2FC| in Table S5")
    unres(claim[0], f"{goid} {desc}", "full GO term membership for a gene-set shift test: resolve from MSigDB "
          "(signature holds only the DEG overlap from Table S7)"
          + ("" if r["Adjusted P-value"] < 0.05 else f"; note term adj. P = {r['Adjusted P-value']:.3g} in Table S7"))

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
print("possible Excel date-corrupted / numeric symbols:", sorted(set(susp.gene) | DATE_GENES))
print("pathway genes missing from DE table:", sig[(sig.kind == "pathway_genes") & (sig.stat == "")].gene.tolist())

claims = yaml.safe_load(open(ROOT / "claims" / f"{PAPER}.yaml"))
need = {c["id"] for c in claims if c["test_template"] in ("state_signature", "geneset_shift", "de_pseudobulk")
        or (c["test_template"].startswith("composition") and (c.get("cell_state") or "").strip())}
in_sig = {x for s in sig.claim_ids.unique() for x in s.split(";")}
in_un = set(un.claim_id)
print(f"\ncoverage: {len(need)} claims need signatures; {len(need & in_sig)} in signatures; "
      f"{len(need & in_un - in_sig)} unresolved only; missing: {sorted(need - in_sig - in_un)}")

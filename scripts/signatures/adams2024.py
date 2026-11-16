"""Signatures for adams2024 (Adams et al. 2024, Nat Aging; SN multiome).

Reads MOESM3_ESM.xlsx (Supplementary Tables 1-6) and writes
claims/signatures/adams2024.tsv and adams2024_unresolved.tsv.
Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/adams2024.py
"""
from pathlib import Path
import pandas as pd
import yaml

PAPER = "adams2024"
ROOT = Path(__file__).resolve().parents[2]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

SUP = cfg.RAW / "supplements" / PAPER
XL = SUP / "43587_2024_583_MOESM3_ESM.xlsx"
OUT = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison",
        "gene", "rank", "stat", "padj", "source", "selection"]
CT = {"N": "Neurons (other)", "ODC": "Oligodendrocytes", "OPC": "OPCs", "AS": "Astrocytes",
      "EC": "Endothelial", "MG": "Microglia", "T": "T cells"}

rows = []


def add(name, claims, kind, ct, state, comp, genes, stats, source, selection):
    genes = [str(g).strip() for g in genes]
    for i, g in enumerate(genes):
        rows.append(dict(signature_id=f"{PAPER}__{name}", paper=PAPER, claim_ids=";".join(claims),
                         kind=kind, cell_type=ct, cell_state=state, comparison=comp, gene=g,
                         rank=i + 1, stat="" if stats is None else stats[i], padj="",
                         source=source, selection=selection))


def col(df, c, start):
    return [g for g in df.iloc[start:, c].tolist() if pd.notna(g) and str(g).strip()]


# ---- Supplementary Table 4: per-gene correlation with cPP score, per glial type --------------
# |r| > 0.1 reproduces the reported counts (ODC 299 up / 474 down; MG 893~894 / 254).
t4 = pd.read_excel(XL, sheet_name="Supplementary Table 4", header=None)
cpp_claims = {"ODC": ["adams2024-C01", "adams2024-C03"], "MG": ["adams2024-C01", "adams2024-C02"],
              "AS": ["adams2024-C01", "adams2024-C17"], "OPC": ["adams2024-C01", "adams2024-C18"]}
for gc in (0, 3, 6, 9):
    ctk = str(t4.iloc[1, gc]).strip()
    d = pd.DataFrame({"gene": t4.iloc[3:, gc], "r": pd.to_numeric(t4.iloc[3:, gc + 1], errors="coerce")}).dropna()
    for kind, sel, extra in [("de_up", d.r > 0.1, ["adams2024-C14"] if ctk == "ODC" else []),
                             ("de_down", d.r < -0.1, ["adams2024-C15"] if ctk == "ODC" else [])]:
        s = d[sel].sort_values("r", ascending=(kind == "de_down"))
        add(f"{ctk}_cPP_{'pos' if kind == 'de_up' else 'neg'}", cpp_claims[ctk] + extra, kind, CT[ctk], "",
            "correlation with cPP score (young -> aged -> PD)", s.gene, s.r.round(4).tolist(),
            f"{XL.name}|Supplementary Table 4",
            f"Pearson r with cPP {'> 0.1' if kind == 'de_up' else '< -0.1'} (reproduces reported counts), ranked by |r|")

# ---- Supplementary Table 3: DEGs per cell type and comparison (gene lists, no stats) ---------
t3 = pd.read_excel(XL, sheet_name="Supplementary Table 3", header=None)
blocks = {0: ("aged_vs_young", "de_down", "aged control vs young control"),   # increased in young vs aged
          8: ("aged_vs_young", "de_up", "aged control vs young control"),
          16: ("PD_vs_aged", "de_up", "PD vs aged control"),
          24: ("PD_vs_aged", "de_down", "PD vs aged control")}
for b0, (cname, kind, comp) in blocks.items():
    assert str(t3.iloc[1, b0]).startswith("Increased in")
    for c in range(b0, b0 + 7):
        ctk = str(t3.iloc[2, c]).strip()
        genes = col(t3, c, 3)
        if not genes:
            continue
        claims = ["adams2024-C26"]
        if ctk == "ODC" and cname == "PD_vs_aged" and kind == "de_down":
            claims += ["adams2024-C24"]
        add(f"{ctk}_{cname}_{kind.split('_')[1]}", claims, kind, CT[ctk], "", comp, genes, None,
            f"{XL.name}|Supplementary Table 3",
            "authors' DEG list (unpaired t-test on nuclei, FC>1.25, P<0.05); order as listed, no stats")

# ---- Supplementary Table 5: ODC state markers and aging/disease gene sets ---------------------
t5 = pd.read_excel(XL, sheet_name="Supplementary Table 5", header=None)
st = {0: ("ODC_healthy_markers", "healthy ODC"), 1: ("ODC_at_risk_markers", "intermediate (at-risk) ODC"),
      2: ("ODC_disease_assoc_markers", "disease-associated ODC (high cPP)")}
for c, (name, state) in st.items():
    add(name, ["adams2024-C04"], "state_markers", "Oligodendrocytes", state, "", col(t5, c, 2), None,
        f"{XL.name}|Supplementary Table 5", "authors' state DEG list (no stats); order as listed, no cap needed (<100)")
assert str(t5.iloc[1, 8]).startswith("Aging Independent") and str(t5.iloc[1, 13]).startswith("Both Aging")
sets = {9: ("ODC_disease_only_gained", "de_up", ["adams2024-C20"]),
        8: ("ODC_disease_only_lost", "de_down", ["adams2024-C19", "adams2024-C24", "adams2024-C25"]),
        14: ("ODC_aging_and_disease_gained", "de_up", ["adams2024-C21"]),
        13: ("ODC_aging_and_disease_lost", "de_down", ["adams2024-C05", "adams2024-C22", "adams2024-C25"])}
for c, (name, kind, claims) in sets.items():
    comp = ("disease-associated vs healthy ODC; aging-independent" if c in (8, 9)
            else "disease-associated vs healthy ODC and aged vs young")
    add(name, claims, kind, "Oligodendrocytes", "disease-associated ODC", comp, col(t5, c, 3), None,
        f"{XL.name}|Supplementary Table 5", "authors' list (overlap of cohort DEGs and ODC-state DEGs); order as listed")

# ---- Named genes from claim text ----------------------------------------------------------------
claims = yaml.safe_load(open(ROOT / "claims" / f"{PAPER}.yaml"))
for c in claims:
    if c["test_template"] in ("de_pseudobulk", "geneset_shift", "state_signature") and c.get("genes"):
        add(f"{c['id'].split('-')[1]}_named", [c["id"]], "named_genes", c["cell_type"], c.get("cell_state") or "",
            c["comparison"], c["genes"], None, "text", "genes named in claim text")

sig = pd.DataFrame(rows, columns=COLS)

# ---- Unresolved ----------------------------------------------------------------------------------
unres = pd.DataFrame([
    ("adams2024-C14", "UPR / chaperone-mediated autophagy module (Methods 'Protein Folding', 'Chaperone-Mediated Autophagy')",
     "module member genes only in Methods (not in supplementary xlsx); resolve centrally from GO/MSigDB. ODC cPP-positive genes attached as proxy"),
    ("adams2024-C15", "myelination / synapse-formation module (Methods 'myelin-forming', 'synaptic support ODC')",
     "module member genes only in Methods; resolve centrally. ODC cPP-negative genes attached as proxy"),
    ("adams2024-C16", "homeostatic MG / stage-1 DAM / aging MG modules (Keren-Shaul 2017 and others)",
     "module gene lists not in supplementary xlsx; take from Keren-Shaul 2017 / Olah 2018 centrally. Only named genes emitted"),
    ("adams2024-C06", "chromatin_da", "no differential peak table in supplement (claim is a null result)"),
    ("adams2024-C07", "chromatin_da", "no promoter accessibility table in supplement"),
    ("adams2024-C08", "chromatin_da peak-gene links (RASGRF1, MSRA)", "peak-gene link tables not in supplement"),
    ("adams2024-C09", "chromatin_da peak-gene links (NEAT1)", "peak-gene link tables not in supplement"),
    ("adams2024-C23", "chromatin_da peak-gene links (MAPT locus)", "peak-gene link tables not in supplement"),
], columns=["claim_id", "term", "reason"])

OUT.mkdir(parents=True, exist_ok=True)
sig.to_csv(OUT / f"{PAPER}.tsv", sep="\t", index=False)
unres.to_csv(OUT / f"{PAPER}_unresolved.tsv", sep="\t", index=False)

# ---- Checks ----------------------------------------------------------------------------------------
assert sig.gene.str.len().gt(0).all()
assert (sig.groupby("signature_id")["rank"].min() == 1).all()
print(sig.groupby(["signature_id", "kind"]).size().to_string())
need = {c["id"] for c in claims if c["test_template"] in ("state_signature", "geneset_shift", "de_pseudobulk")
        or (c["test_template"] == "composition_donor" and c.get("cell_state"))}
chrom = {c["id"] for c in claims if c["test_template"] == "chromatin_da"}
covered = set(";".join(sig.claim_ids).split(";"))
unr = set(unres.claim_id)
print(f"\nrequired claims: {len(need)}; in signatures: {len(need & covered)}; unresolved only: {len(need & unr - covered)}")
print("MISSING:", sorted(need - covered - unr) or "none")
print("chromatin claims unresolved:", sorted(chrom & unr), "missing:", sorted(chrom - unr - covered) or "none")

"""Signatures for feleke2021 (Feleke et al. 2021, Acta Neuropathol; ACC snRNA, 7 CTRL / 7 PD / 7 PDD / 7 DLB).

Reads (Supplementary Table numbering = MOESM number - 1):
  401_2021_2343_MOESM6_ESM.xlsx  Supplementary Table 5: snRNA DEGs (MAST), all pairwise comparisons
  401_2021_2343_MOESM7_ESM.xlsx  Supplementary Table 6: GO enrichment of DEG sets (overlapping genes)
  401_2021_2343_MOESM8_ESM.xlsx  Supplementary Table 7: PD PRS-derived pathway (Bandres-Ciga 2020) enrichment
and writes claims/signatures/feleke2021.tsv and claims/signatures/feleke2021_unresolved.tsv.
Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/feleke2021.py

Notes
- Comparison "A_vs_B": A is the reference; coef = log2FC of B relative to A (column description).
  Checked: SNCA coef > 0 in Control_vs_PD for ExN/InN/astro/oligo, as in the text (claim C11).
- The ci_high / ci_low columns are inconsistent with coef (sign opposite to coef in every row, coef
  inside the interval in only ~12% of rows); they are ignored. Direction is taken from coef, which
  matches the authors' up/down labels in Supplementary Table 6.
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import yaml

PAPER = "feleke2021"
ROOT = Path(__file__).resolve().parents[2]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

SUP = cfg.RAW / "supplements" / PAPER
X_DE = SUP / "401_2021_2343_MOESM6_ESM.xlsx"
X_GO = SUP / "401_2021_2343_MOESM7_ESM.xlsx"
X_PRS = SUP / "401_2021_2343_MOESM8_ESM.xlsx"
OUT = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison",
        "gene", "rank", "stat", "padj", "source", "selection"]
C = lambda *n: [f"{PAPER}-C{int(i):02d}" for i in n]  # noqa: E731
CT = {"Astrocyte": ("Astrocytes", "Astro"), "Excitatory neuron": ("Excitatory neurons", "ExN"),
      "Inhibitory neuron": ("Inhibitory neurons", "InN"), "Microglia": ("Microglia", "Micro"),
      "OPC": ("OPCs", "OPC"), "Oligodendrocyte": ("Oligodendrocytes", "Oligo"),
      "Vascular": ("Endothelial", "Vasc")}
LBL = {"Control": "control", "PD": "PD", "PDD": "PDD", "DLB": "DLB"}

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


# ---------------- Supplementary Table 5: DEGs ----------------
de = pd.read_excel(X_DE, sheet_name="Results")
n0 = len(de)
de = de[(de.fdr < 0.05) & (de.coef.abs() > np.log2(1.5))]
print(f"ST5 rows {n0}; after fdr<0.05 & |log2FC|>log2(1.5): {len(de)}")
claims_by_comp = {"Control_vs_PD": C(6, 7, 10, 13), "Control_vs_PDD": C(8, 13), "Control_vs_DLB": C(8, 13),
                  "PD_vs_PDD": C(9), "PD_vs_DLB": C(9), "PDD_vs_DLB": C(9)}
for (comp, ct), d in de.groupby(["comparison", "cell_type"]):
    ref, grp = comp.split("_vs_")
    comp_lbl = f"{LBL[grp]} vs {LBL[ref]}"
    ctn, cts = CT[ct]
    cl = list(claims_by_comp[comp])
    if comp == "Control_vs_PD" and ct in ("Excitatory neuron", "Inhibitory neuron", "Astrocyte", "Oligodendrocyte"):
        cl.append(C(11)[0])
    if comp in ("Control_vs_PDD", "Control_vs_DLB") and ct == "Excitatory neuron":
        cl.append(C(12)[0])
    state = "vascular cells (endothelial + pericytes)" if ct == "Vascular" else ""
    for kind, sub in (("de_up", d[d.coef > 0]), ("de_down", d[d.coef < 0])):
        sub = sub.reindex(sub.coef.abs().sort_values(ascending=False).index)
        if len(sub):
            add(f"{cts}_{grp}_vs_{ref}_{kind.split('_')[1]}", sorted(set(cl)), kind, ctn, state, comp_lbl,
                sub.hgnc_symbol, sub.coef.round(4).tolist(), sub.fdr.tolist(), f"{X_DE.name}|Results",
                "Supplementary Table 5 (MAST, cells): fdr<0.05 & |coef (log2FC)|>log2(1.5); ranked by |log2FC|")

# ---------------- Supplementary Table 6: GO terms named in claim C10 ----------------
go = pd.read_excel(X_GO, sheet_name="Results")
go_pd = go[go.comparison == "Control_vs_PD"]
for term, short in (("GO:0098978", "glutamatergic_synapse"), ("GO:0005743", "mito_inner_membrane"),
                    ("GO:0043687", "post_translational_protein_modification")):
    t = go_pd[go_pd.parent_term == term]
    desc = t.parent_description.iloc[0]
    genes = pd.Series([g for s in t.overlapping_genes.dropna() for g in str(s).split(";") if g.strip()])
    cnt = genes.value_counts()
    cnt = cnt.reindex(sorted(cnt.index, key=lambda g: (-cnt[g], g)))
    add(f"GO_{short}_PD_vs_control", C(10), "pathway_genes", "All", "", "PD vs control",
        cnt.index, cnt.tolist(), None, f"{X_GO.name}|Results",
        f"union of overlapping DEGs of all child terms under parent {term} ({desc}), PD vs control, all cell "
        f"types and both directions ({t.cell_type.nunique()} cell types); stat = number of cell type x direction x "
        "child-term sets containing the gene")
    unres(C(10)[0], f"{term} {desc}", "full GO term membership for a gene-set shift test: resolve from MSigDB "
          "(signature holds only the DEG overlap reported in Supplementary Table 6)")

# ---------------- named genes ----------------
add("SNCA_PD_up", C(11), "named_genes", "All", "ExN, InN, astrocytes, oligodendrocytes", "PD vs control",
    ["SNCA"], None, None, "text", "text list")
add("SNCA_ExN", C(12), "named_genes", "Excitatory neurons", "SNCA-high excitatory neurons", "PDD vs control; DLB vs control",
    ["SNCA"], None, None, "text", "text list")

# ---------------- unresolved ----------------
unres(C(5)[0], "vascular cells (endothelial + pericytes merged)",
      "composition of a broad merged class (no marker table in supplement); use endothelial + pericyte annotation; "
      "proportions in Supplementary Table 2 are Scaden deconvolution of bulk RNA-seq")
prs = pd.read_excel(X_PRS, sheet_name="Results").drop_duplicates("pathway")
for _, r in prs.iterrows():
    link = str(r.link)
    pid = link.rsplit("/", 1)[-1] if "msigdb" in link else link
    unres(C(13)[0], f"{pid} {r.pathway} ({r.pathway_group})",
          "PRS-derived pathway (Bandres-Ciga 2020); full membership not in supplement (only DEG overlaps); "
          "resolve from Reactome/KEGG via MSigDB; DEG lists vs control attached")

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
susp = sig[sig.gene.str.match(r"^\d+-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$|^(Jan|Feb|Mar|Apr|Sep|Oct|Nov|Dec)-\d+$")]
print("possible Excel date-corrupted symbols:", sorted(set(susp.gene) | DATE_GENES))

claims = yaml.safe_load(open(ROOT / "claims" / f"{PAPER}.yaml"))
need = {c["id"] for c in claims if c["test_template"] in ("state_signature", "geneset_shift", "de_pseudobulk")
        or (c["test_template"].startswith("composition") and (c.get("cell_state") or "").strip())}
in_sig = {x for s in sig.claim_ids.unique() for x in s.split(";")}
in_un = set(un.claim_id)
print(f"\ncoverage: {len(need)} claims need signatures; {len(need & in_sig)} in signatures; "
      f"{len(need & in_un - in_sig)} unresolved only; missing: {sorted(need - in_sig - in_un)}")

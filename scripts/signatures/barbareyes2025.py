"""Gene signatures for barbareyes2025 (oligodendroglia in the PD dorsal striatum, Acta Neuropathol 2025).

Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/barbareyes2025.py
Writes claims/signatures/barbareyes2025.tsv and claims/signatures/barbareyes2025_unresolved.tsv (claims/SIGNATURES.md).

Source: Online Resource 2 (401_2025_2884_MOESM2_ESM.xlsx).
- ST2 '<state>_markers': scanpy Wilcoxon markers per subclass. Each gene appears once, with the comparison column
  '<state> vs:' naming one other subclass (pairwise tests collapsed by the authors), so the markers mix lineage-level
  genes (e.g. vs OPCs) with state-specific ones. The text marker genes are emitted separately as `*_markers_text`.
- ST3 GSEA (GSEApy prerank, cluster vs cluster): leading-edge genes of significant terms, grouped by the themes named in
  the claims.
"""
from pathlib import Path
import re
import numpy as np
import pandas as pd
import yaml

PAPER = "barbareyes2025"
ROOT = Path(__file__).resolve().parents[2]
OUTDIR = ROOT / "claims" / "signatures"

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

SUP = cfg.RAW / "supplements" / PAPER
XLSX = SUP / "401_2025_2884_MOESM2_ESM.xlsx"
SRC = XLSX.name
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison", "gene", "rank", "stat",
        "padj", "source", "selection"]
rows, unres = [], []


def cid(*nums):
    return [f"{PAPER}-C{n:02d}" for n in nums]


def add(name, claims, kind, cell_type, state, comparison, genes, source, selection, stat=None, padj=None):
    sid = f"{PAPER}__{name}"
    seen, n = set(), 0
    for i, g in enumerate(genes):
        if g is None or (isinstance(g, float) and np.isnan(g)) or str(g).strip() == "" or g in seen:
            continue
        seen.add(g); n += 1
        rows.append(dict(signature_id=sid, paper=PAPER, claim_ids=";".join(claims), kind=kind, cell_type=cell_type,
                         cell_state=state, comparison=comparison, gene=str(g).strip(), rank=n,
                         stat="" if stat is None else stat[i], padj="" if padj is None else padj[i],
                         source=source, selection=selection))
    return n


def unresolved(claims, term, reason):
    for c in claims:
        unres.append(dict(claim_id=c, term=term, reason=reason))


x = pd.ExcelFile(XLSX)

# ---------------------------------------------------------------- state markers (ST2)
STATES = {  # state: (sheet, cell_type, claims, text markers)
    "PDAOPCs": ("ST2. PDAOPCs_markers", "OPCs", cid(1, 6), ["VCAN", "MT3", "CD81"]),
    "OPCs-A": ("ST2. OPCs-A_markers", "OPCs", [], []),
    "PDAO-1": ("ST2. PDAO1_markers", "Oligodendrocytes", cid(1, 3, 10), ["OPALIN", "DBNDD2", "FTL", "HSPA1A"]),
    "PDAO-2": ("ST2. PDAO2_markers", "Oligodendrocytes", cid(1, 4, 11), ["CTNNA2", "ARHGAP24"]),
    "PDAO-3": ("ST2. PDAO3_markers", "Oligodendrocytes", cid(1, 5), ["SLC5A11", "ACTN2", "MT-ND3"]),
    "MOL-A": ("ST2. MOL-A_markers", "Oligodendrocytes", cid(7), []),
    "MOL-F": ("ST2. MOL-F_markers", "Oligodendrocytes", cid(8), []),
    "MOL-G": ("ST2. MOL-G_markers", "Oligodendrocytes", [], []),
}
for st, (sheet, ct, cl, text) in STATES.items():
    d = x.parse(sheet)
    s = d[(d.pvals_adj < 0.05) & (d.logfoldchanges > 0)].sort_values("logfoldchanges", ascending=False)
    s = s.drop_duplicates("names").head(100)
    nm = st.replace("-", "")
    add(f"{nm}_markers", cl, "state_markers", ct, st, "", s.names.tolist(), f"{SRC}|{sheet}",
        f"padj<0.05 & logFC>0, top 100 by logFC (Wilcoxon vs other subclasses; table keeps one pairwise comparison "
        f"per gene, so lineage-level genes are included)", stat=s.logfoldchanges.round(4).tolist(),
        padj=s.pvals_adj.tolist())
    if text:
        add(f"{nm}_markers_text", cl, "state_markers", ct, st, "", text, "text", "text list (markers named in Results)")

# ---------------------------------------------------------------- named genes
NAMED = [
    # DNAJB1 in the claim's gene list was added by the curator ('DNAJ family members' in the text) and is not used here
    ("HSP_genes_PDA_states", cid(2), "Oligodendrocytes", "PDAOPCs; PDAO-1; PDAO-2; PDAO-3",
     "PD-associated subpopulations vs other oligodendroglia subpopulations",
     ["HSPA1A", "HSPB1", "HSP90AB1", "HSPA1B", "HSPA8", "HSP90AA1", "HSPA5", "CRYAB"]),
    ("PDAO1_hotspot_genes", cid(10), "Oligodendrocytes", "PDAO-1", "PDAO-1 vs other OPALIN+ MOLs",
     ["CRYAB", "FTH1", "FTL", "HSPA1A", "HSPA1B", "HSPH1"]),
    ("PDAO2_hotspot_genes", cid(11), "Oligodendrocytes", "PDAO-2", "PDAO-2 vs other OPALIN+ MOLs",
     ["ETV5", "QDPR", "ARHGAP24", "CTNNA2"]),
    ("GABRR2_PDAO1_PDAO2", cid(25), "Oligodendrocytes", "PDAO-1; PDAO-2",
     "PDAO-1/PDAO-2 vs other oligodendroglia subpopulations", ["GABRR2"]),
    ("PDAO2_receptor_genes", cid(26), "Oligodendrocytes", "PDAO-2", "PDAO-2 vs other oligodendroglia subpopulations",
     ["CHRNA10", "GABRG1"]),
    ("HSPR_module_xenium", cid(29, 31, 32), "Oligodendrocytes", "HSPR module", "PD vs control (Xenium)",
     ["HSPA1A", "DNAJA1", "HSP90AB1", "CRYAB"]),
    ("myelination_module_xenium", cid(30, 31), "Oligodendrocytes", "myelination module", "PD vs control (Xenium)",
     ["MAG", "MOBP", "CNTN2"]),
]
for name, cl, ct, st, comp, genes in NAMED:
    add(name, cl, "named_genes", ct, st, comp, genes, "text", "genes named in the text (order as written)")

# ---------------------------------------------------------------- GSEA leading edges (ST3)
THEMES = {
    "translation_ribosome": r"\bTranslation \(|Cytoplasmic Translation|Ribosom|rRNA",
    "protein_folding_stress": r"Protein Folding|Refolding|Unfolded Protein|Heat|Chaperone|Endoplasmic Reticulum Stress|Integrated Stress",
    "mitochondrial_respiration": r"Mitochondri|Oxidative Phosphorylation|Respiration|Electron Transport",
    "apoptosis": r"Apopto|Cell Death",
    "immune": r"Immun|Interferon|Antigen|Cytokine",
}
COMPS = [("ST3. PDAO-1 vs MOL-A", "PDAO-1", "MOL-A", "Oligodendrocytes"),
         ("ST3. PDAO-2 vs MOL-A", "PDAO-2", "MOL-A", "Oligodendrocytes"),
         ("ST3. PDAO-3 vs MOL-F", "PDAO-3", "MOL-F", "Oligodendrocytes"),
         ("ST3. PDAOPCs vs OPCs-A", "PDAOPCs", "OPCs-A", "OPCs")]


def leading_edge(d, pattern):
    m = d[d.Term.str.contains(pattern, case=False, regex=True)]
    cnt = {}
    for _, r in m.iterrows():
        for g in str(r.Lead_genes).split(";"):
            g = g.strip()
            if g:
                cnt[g] = cnt.get(g, 0) + 1
    s = pd.Series(cnt, dtype=float).sort_values(ascending=False, kind="stable")
    return s, len(m)


GSEA = {}
for sheet, pdst, ref, ct in COMPS:
    d = x.parse(sheet)
    d = d[(d["FDR q-val"] < 0.05) & (d["UP-regulated in:"] == pdst)]
    GSEA[pdst] = d
    for theme, pat in THEMES.items():
        s, nterm = leading_edge(d, pat)
        nm = f"{pdst.replace('-', '')}_vs_{ref.replace('-', '')}_{theme}"
        n = add(nm, cid(12), "pathway_genes", ct, pdst, f"{pdst} vs {ref}", s.index.tolist(), f"{SRC}|{sheet}",
                f"leading-edge genes of GSEA terms with FDR q<0.05 up in {pdst} whose name matches /{pat}/ "
                f"({nterm} terms); ranked by number of matching terms (stat)", stat=s.astype(int).tolist())
        if n == 0:
            unresolved(cid(12), nm, f"no significant {theme} term up in {pdst} vs {ref}")

# C13: PDAO-2 DAP12 / interferon
s, nterm = leading_edge(GSEA["PDAO-2"], r"DAP12|Interferon")
add("PDAO2_vs_MOLA_DAP12_interferon", cid(13), "pathway_genes", "Oligodendrocytes", "PDAO-2", "PDAO-2 vs MOL-A",
    s.index.tolist(), f"{SRC}|ST3. PDAO-2 vs MOL-A", f"leading-edge genes of GSEA terms with FDR q<0.05 up in PDAO-2 "
    f"matching /DAP12|Interferon/ ({nterm} terms); ranked by number of matching terms", stat=s.astype(int).tolist())
# C14: PDAOPCs adaptive immune, antigen receptor, DNA repair
for theme, pat in [("adaptive_immune_antigen_receptor", r"T Cell|B Cell|Lymphocyte|Adaptive|Antigen"),
                   ("DNA_repair", r"DNA Repair|Double-Strand Break|Excision Repair|Mismatch Repair")]:
    s, nterm = leading_edge(GSEA["PDAOPCs"], pat)
    nm = f"PDAOPCs_vs_OPCsA_{theme}"
    n = add(nm, cid(14), "pathway_genes", "OPCs", "PDAOPCs", "PDAOPCs vs OPCs-A", s.index.tolist(),
            f"{SRC}|ST3. PDAOPCs vs OPCs-A", f"leading-edge genes of GSEA terms with FDR q<0.05 up in PDAOPCs matching "
            f"/{pat}/ ({nterm} terms); ranked by number of matching terms", stat=s.astype(int).tolist())
    if n == 0:
        unresolved(cid(14), nm, "no significant matching term")

# ---------------------------------------------------------------- unresolved
unresolved(cid(9), "all 15 oligodendroglia subpopulations",
           "composition null claim (sex independence) over all 15 subclasses; needs label transfer of every subclass; "
           "marker signatures emitted only for PDAOPCs, OPCs-A, PDAO-1/2/3, MOL-A/F/G")
unresolved(cid(11), "GO: response to cytokine; antigen processing and presentation (Hotspot modules 6, 13)",
           "module gene lists only in Supplementary Fig. 5 (not tabulated); named genes used")
unresolved(cid(10), "Hotspot modules 1 and 2 (PDAO-1)", "full module lists only in Supplementary Fig. 5; named genes used")
unresolved(cid(29, 30, 32), "Xenium Hotspot modules 2 (HSPR) and 3 (myelination)",
           "full module lists not in the supplement; genes named in Fig. 6d used")

# ---------------------------------------------------------------- write + checks
sig = pd.DataFrame(rows, columns=COLS)
un = pd.DataFrame(unres, columns=["claim_id", "term", "reason"])
OUTDIR.mkdir(parents=True, exist_ok=True)
sig.to_csv(OUTDIR / f"{PAPER}.tsv", sep="\t", index=False)
un.to_csv(OUTDIR / f"{PAPER}_unresolved.tsv", sep="\t", index=False)

assert sig.gene.str.strip().ne("").all() and sig.gene.notna().all()
assert (sig.groupby("signature_id")["rank"].min() == 1).all()
assert not sig.duplicated(["signature_id", "gene"]).any()
summ = sig.groupby(["signature_id", "kind"]).size().rename("n_genes").reset_index()
pd.set_option("display.width", 200); pd.set_option("display.max_rows", 500)
print(summ.to_string(index=False))
print(f"\n{sig.signature_id.nunique()} signatures, {len(sig)} rows; by kind:")
print(summ.kind.value_counts().to_string())

claims = yaml.safe_load(open(ROOT / "claims" / f"{PAPER}.yaml"))
need = [c["id"] for c in claims if c["test_template"] in ("state_signature", "geneset_shift", "de_pseudobulk")
        or ((c["test_template"] == "composition_donor" or c["claim_type"] == "composition") and (c.get("cell_state") or "").strip())]
in_sig = set(";".join(sig.claim_ids).split(";"))
in_un = set(un.claim_id)
missing = [c for c in need if c not in in_sig and c not in in_un]
print(f"\nCoverage: {len(need)} claims need a signature; {sum(c in in_sig for c in need)} in signatures, "
      f"{sum(c not in in_sig and c in in_un for c in need)} only in unresolved, missing: {missing or 'none'}")

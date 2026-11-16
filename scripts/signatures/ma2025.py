"""Signatures for ma2025 (Ma et al. 2025, Nat Commun; SN + cingulate snRNA, Visium, TCR).

MOESM numbering: Supplementary Data n = MOESM(n+2) (from MOESM2 description file).
  SD4 = MOESM6 (SN DEGs, pathfindR, cluster markers)   SD5 = MOESM7 (cingulate DEGs, pathfindR)
  SD6 = MOESM8 (scHPF factors; sheet C = T-cell gene scores)   SD7 = MOESM9 (TRM / memory gene sets; strict OOXML)
  SD9 = MOESM11 (Visium DEGs: A = SNpc, B = surrounding tissue)
Writes claims/signatures/ma2025.tsv and ma2025_unresolved.tsv.
Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/ma2025.py
"""
from pathlib import Path
import re
import warnings
import zipfile
import pandas as pd
import yaml

warnings.filterwarnings("ignore")
PAPER = "ma2025"
ROOT = Path(__file__).resolve().parents[2]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

SUP = cfg.RAW / "supplements" / PAPER
F = lambda m: SUP / f"41467_2025_62478_MOESM{m}_ESM.xlsx"
OUT = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison",
        "gene", "rank", "stat", "padj", "source", "selection"]
rows = []


def add(name, claims, kind, ct, state, comp, genes, stats, padj, source, selection):
    if stats is None and padj is None:
        genes = list(dict.fromkeys(str(g).strip() for g in genes))  # de-duplicate text lists (SD7 B lists ITGB1 twice)
    for i, g in enumerate(genes):
        rows.append(dict(signature_id=f"{PAPER}__{name}", paper=PAPER, claim_ids=";".join(claims), kind=kind,
                         cell_type=ct, cell_state=state, comparison=comp, gene=str(g).strip(), rank=i + 1,
                         stat="" if stats is None else stats[i], padj="" if padj is None else padj[i],
                         source=source, selection=selection))


def de_blocks(m, sheet):
    """Side-by-side limma blocks: [gene, logFC, AveExpr, t, P.Value, adj.P.Val, B, new.adj.P.Val]."""
    d = pd.read_excel(F(m), sheet_name=sheet, header=None)
    hr = 1 if "logFC" in d.iloc[1].astype(str).tolist() else 0
    h = d.iloc[hr].astype(str).tolist()
    out = {}
    for st in [i - 1 for i, v in enumerate(h) if v == "logFC"]:
        name = "NA"
        if hr == 1:
            for c in (st, st + 1):
                if pd.notna(d.iloc[0, c]):
                    name = str(d.iloc[0, c]).strip()
                    break
        b = d.iloc[hr + 1:, st:st + 8].copy()
        b.columns = ["gene", "logFC", "AveExpr", "t", "P", "adjP", "B", "newadjP"]
        b = b.dropna(subset=["gene"])
        for c in b.columns[1:]:
            b[c] = pd.to_numeric(b[c], errors="coerce")
        out[name] = b
    return out


def path_blocks(m, sheet):
    d = pd.read_excel(F(m), sheet_name=sheet, header=None)
    hr = [i for i in range(3) if "Term_Description" in d.iloc[i].astype(str).tolist()][0]
    h = d.iloc[hr].astype(str).tolist()
    out = {}
    for st in [i - 1 for i, v in enumerate(h) if v == "ID"]:
        name = "NA"
        if hr >= 1:
            for c in (st, st + 1, st - 1):
                if 0 <= c < d.shape[1] and pd.notna(d.iloc[0, c]):
                    name = str(d.iloc[0, c]).strip()
                    break
        b = d.iloc[hr + 1:, st + 1:st + 10].copy()
        b.columns = ["ID", "term", "FE", "occ", "support", "lowest_p", "highest_p", "up", "down"]
        out[name] = b.dropna(subset=["ID"])
    return out


def strict_xlsx_col_a(path):
    """MOESM9 is strict-OOXML (openpyxl sees no sheets): read column A of each sheet from the XML."""
    z = zipfile.ZipFile(path)
    ss = [re.sub(r"<[^>]+>", "", s) for s in re.findall(r"<si>(.*?)</si>", z.read("xl/sharedStrings.xml").decode(), re.S)]
    names = re.findall(r'<sheet name="([^"]+)"', z.read("xl/workbook.xml").decode())
    out = {}
    for i, n in enumerate(names, 1):
        x = z.read(f"xl/worksheets/sheet{i}.xml").decode()
        cells = re.findall(r'<c r="A(\d+)"[^>]*t="s"[^>]*><v>(\d+)</v>', x)
        vals = [ss[int(v)].strip() for r, v in sorted(cells, key=lambda t: int(t[0]))]
        out[n] = [v for v in vals[1:] if v]  # drop header 'Gene'
    return out


DE_SEL = "authors' limma DE (cells; donor/age/sex covariates): adj.P.Val<0.05, split by sign of logFC, ranked by logFC"
CT = {"Protoplasmic": ("Astrocytes", "protoplasmic astrocytes"), "Fibrous": ("Astrocytes", "fibrous-like astrocytes"),
      "Dopaminergic": ("DA neurons", ""), "MEG3 Interneurons": ("Neurons (other)", "MEG3 interneurons"),
      "GABAergic": ("Inhibitory neurons", ""), "Quiescent": ("Microglia", "quiescent microglia"),
      "Activated": ("Microglia", "activated microglia"), "BAM": ("Immune (other)", "BAMs"),
      "Oligo": ("Oligodendrocytes", ""), "OPC": ("OPCs", ""), "Tcells": ("T cells", ""),
      "Vascular": ("Pericytes", "vascular (mural)"), "Endothelial": ("Endothelial", ""),
      "VascEndo": ("Endothelial", "vascular and endothelial (pooled)")}

# claims attached to DE lists: (region, block, direction) -> claims
DE_CLAIMS = {("SN", "Dopaminergic", "down"): ["ma2025-C09", "ma2025-C10"], ("SN", "Dopaminergic", "up"): ["ma2025-C10"],
             ("SN", "Activated", "up"): ["ma2025-C11", "ma2025-C13", "ma2025-C15"],
             ("SN", "Activated", "down"): ["ma2025-C12", "ma2025-C14", "ma2025-C15"],
             ("CC", "Activated", "up"): ["ma2025-C11", "ma2025-C15"], ("CC", "Activated", "down"): ["ma2025-C12", "ma2025-C15"],
             ("SN", "Protoplasmic", "up"): ["ma2025-C18", "ma2025-C23", "ma2025-C25", "ma2025-C27"],
             ("SN", "Protoplasmic", "down"): ["ma2025-C20"],
             ("CC", "Protoplasmic", "up"): ["ma2025-C19", "ma2025-C21"], ("CC", "Protoplasmic", "down"): ["ma2025-C22"],
             ("SN", "Tcells", "up"): ["ma2025-C08"]}

for region, m, sheets in [("SN", 6, "ACEGHJ"), ("CC", 7, "ACEGH")]:
    comp = "PD vs control (SN)" if region == "SN" else "PDD/DLBD vs control (cingulate)"
    for s in sheets:
        for name, b in de_blocks(m, s).items():
            if name == "NA":
                name = "Tcells" if (m, s) == (6, "H") else "VascEndo"  # per MOESM2 sheet descriptions
            if name in CT:
                ct, state = CT[name]
            elif name.endswith("_Glu"):
                ct, state = "Excitatory neurons", name
            elif name.endswith("_GABA"):
                ct, state = "Inhibitory neurons", name
            else:
                raise ValueError(name)
            sig = b[b.adjP < 0.05]
            for d, sub in (("up", sig[sig.logFC > 0].sort_values("logFC", ascending=False)),
                           ("down", sig[sig.logFC < 0].sort_values("logFC"))):
                if sub.empty:
                    continue
                add(f"{region}_{name.replace(' ', '_')}_{d}", DE_CLAIMS.get((region, name, d), []), f"de_{d}", ct, state,
                    comp, sub.gene.tolist(), sub.logFC.round(4).tolist(), sub.adjP.tolist(),
                    f"{F(m).name}|{s}", DE_SEL)

# ---- SN cluster markers (SD4 sheet K, Seurat FindAllMarkers) ---------------------------------------------
k = pd.read_excel(F(6), sheet_name="K")
mk_claims = {"Myeloid_Activated": (["ma2025-C16"], "Microglia", "CD68+ activated microglia"),
             "Fibrous": (["ma2025-C23", "ma2025-C24"], "Astrocytes", "CD44+ fibrous-like astrocytes"),
             "Myeloid_BAM": (["ma2025-C17"], "Immune (other)", "BAMs"),
             "Tcells": (["ma2025-C01"], "T cells", "T cells (CD8+ dominated)")}
for cl, (cids, ct, state) in mk_claims.items():
    g = k[(k.cluster == cl) & (k.p_val_adj < 0.05) & (k.avg_log2FC > 0)].sort_values("avg_log2FC", ascending=False).head(100)
    add(f"SN_{cl}_markers", cids, "state_markers", ct, state, "", g.gene.tolist(), g.avg_log2FC.round(4).tolist(),
        g.p_val_adj.tolist(), f"{F(6).name}|K", "p_val_adj<0.05 & avg_log2FC>0, top 100 by avg_log2FC (cluster vs rest, SN)")

# ---- T-cell scHPF factors (SD6 sheet C) ----------------------------------------------------------------
tc = pd.read_excel(F(8), sheet_name="C")
for fac, cids, state in [("Factor5", ["ma2025-C06", "ma2025-C26"], "T cell Factor 5 (TRM)"),
                         ("Factor1", ["ma2025-C07"], "T cell Factor 1 (effector)")]:
    g = tc.sort_values(fac, ascending=False).head(100)
    add(f"Tcell_scHPF_{fac}", cids, "state_markers", "T cells", state, "", g.Gene.tolist(), g[fac].round(5).tolist(),
        None, f"{F(8).name}|C", f"top 100 genes by scHPF gene score ({fac}); no p-values")

# ---- TRM / memory gene sets (SD7) ------------------------------------------------------------------------
gs = strict_xlsx_col_a(F(9))
add("CD8_TRM_geneset", ["ma2025-C06", "ma2025-C08", "ma2025-C26"], "state_markers", "T cells", "CD8+ TRM", "",
    gs["A"], None, None, f"{F(9).name}|A", "text list (authors' curated CD8+ TRM gene set; order as listed)")
add("general_memory_geneset", ["ma2025-C06"], "state_markers", "T cells", "memory T cells", "",
    gs["B"], None, None, f"{F(9).name}|B", "text list (authors' curated general memory gene set; order as listed)")

# ---- Visium DE (SD9) ------------------------------------------------------------------------------------
for s, area, cids in [("A", "SNpc", {"up": ["ma2025-C29"], "down": ["ma2025-C28", "ma2025-C29"]}),
                      ("B", "surrounding", {"up": ["ma2025-C29", "ma2025-C30"], "down": ["ma2025-C29"]})]:
    v = pd.read_excel(F(11), sheet_name=s)
    v = v[v.p_val_adj < 0.05]
    for d, sub in (("up", v[v.log2FC > 0].sort_values("log2FC", ascending=False)),
                   ("down", v[v.log2FC < 0].sort_values("log2FC"))):
        add(f"Visium_{area}_{d}", cids[d], f"de_{d}", "All", f"Visium {area} spots", f"PD vs control (Visium {area})",
            sub.Gene.tolist(), sub.log2FC.round(4).tolist(), sub.p_val_adj.tolist(), f"{F(11).name}|{s}",
            "p_val_adj<0.05 (limma on spots), split by sign, ranked by log2FC")

# ---- pathfindR pathway member genes ---------------------------------------------------------------------
PW = [  # (claim, MOESM, sheet, block, KEGG ids, cell_type, state, region)
    ("ma2025-C10", 6, "D", "Dopaminergic", ["hsa00350", "hsa00410"], "DA neurons", "", "SN"),
    ("ma2025-C15", 6, "F", "Activated", ["hsa04612", "hsa04659", "hsa04658", "hsa04660", "hsa04623"], "Microglia", "activated microglia", "SN"),
    ("ma2025-C15", 7, "F", "Activated", ["hsa04612", "hsa04659", "hsa04658", "hsa04660"], "Microglia", "activated microglia", "CC"),
    ("ma2025-C25", 6, "B", "Protoplasmic", ["hsa04630", "hsa04064", "hsa04657", "hsa04010"], "Astrocytes", "protoplasmic astrocytes", "SN"),
]
pw_found = set()
for cid, m, s, blk, ids, ct, state, region in PW:
    b = path_blocks(m, s)[blk]
    for kid in ids:
        r = b[b.ID == kid]
        if r.empty:
            continue
        pw_found.add((cid, kid))
        r = r.iloc[0]
        split = lambda x: [g.strip() for g in str(x).split(",") if g.strip() and str(x).strip() != "nan"]
        up, dn = split(r.up), split(r.down)
        add(f"{region}_{blk}_{kid}", [cid], "pathway_genes", ct, state,
            "PD vs control (SN)" if region == "SN" else "PDD/DLBD vs control (cingulate)",
            up + dn, ["up"] * len(up) + ["down"] * len(dn), None, f"{F(m).name}|{s}",
            f"pathfindR {kid} {r.term} (lowest_p={float(r.lowest_p):.2g}); members as listed, up first; stat = DE direction")

# ---- named genes -----------------------------------------------------------------------------------------
claims = yaml.safe_load(open(ROOT / "claims" / f"{PAPER}.yaml"))
for c in claims:
    if (c["test_template"] in ("de_pseudobulk", "state_signature", "geneset_shift")
            or (c["test_template"] == "composition_donor" and c.get("cell_state"))) and c.get("genes"):
        add(f"{c['id'].split('-')[1]}_named", [c["id"]], "named_genes", c["cell_type"], c.get("cell_state") or "",
            c["comparison"], c["genes"], None, None, "text", "genes named in claim text")

sig = pd.DataFrame(rows, columns=COLS)

unres = [("ma2025-C17", "Schonhoff et al. 2023 murine proliferative/activated BAM signatures",
          "external mouse gene sets (ref. 24), not in supplement; SN BAM cluster markers attached; SN BAM DE has 0 genes at adj.P<0.05"),
         ("ma2025-C24", "astrocyte gene programs 8-10 (Fig. S10)",
          "program gene loadings not in supplement (SD6 sheet A factors 8-10 are not astrocytic); SN fibrous-like astrocyte markers attached")]
for cid, m, s, blk, ids, *_ in PW:
    for kid in ids:
        if (cid, kid) not in pw_found:
            unres.append((cid, f"KEGG {kid}", f"term not in authors' pathfindR table MOESM{m}|{s} block {blk}; resolve from MSigDB"))
unres = pd.DataFrame(unres, columns=["claim_id", "term", "reason"])

OUT.mkdir(parents=True, exist_ok=True)
sig.to_csv(OUT / f"{PAPER}.tsv", sep="\t", index=False)
unres.to_csv(OUT / f"{PAPER}_unresolved.tsv", sep="\t", index=False)

assert sig.gene.str.len().gt(0).all() and (sig.groupby("signature_id")["rank"].min() == 1).all()
print(sig.groupby(["signature_id", "kind"]).size().to_string())
print(unres.to_string())
need = {c["id"] for c in claims if c["test_template"] in ("state_signature", "geneset_shift", "de_pseudobulk")
        or (c["test_template"] == "composition_donor" and c.get("cell_state"))}
covered = set(";".join(sig.claim_ids).split(";"))
print(f"\nrequired claims: {len(need)}; in signatures: {len(need & covered)}; MISSING:",
      sorted(need - covered - set(unres.claim_id)) or "none")

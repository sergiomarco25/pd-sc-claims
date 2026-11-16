"""Gene signatures for lee2026 (PsychAD cross-disorder atlas, DLPFC; Nature).

Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/lee2026.py
Writes claims/signatures/lee2026.tsv and claims/signatures/lee2026_unresolved.tsv (claims/SIGNATURES.md).

Sources: 41586_2025_9573_MOESM3_ESM.zip (Supplementary Tables), extracted to MOESM3_ESM/ in the same folder.
- Table 4: dreamlet pseudobulk DE per subclass and disorder (sheets PD, DLBD, ...; columns estimate = logFC, FDR,
  crossDis_PP / crossDis_cat = shared cross-disorder component) and sheet 'Cross-disorder Shared Genes' (no header).
- Table 10: Xenium custom panel markers per class/subclass/subtype (the only per-subclass marker list supplied;
  Table 2 is the taxonomy without genes).
"""
from pathlib import Path
import zipfile
import numpy as np
import pandas as pd
import yaml

PAPER = "lee2026"
ROOT = Path(__file__).resolve().parents[2]
OUTDIR = ROOT / "claims" / "signatures"

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

SUP = cfg.RAW / "supplements" / PAPER
ZIP = SUP / "41586_2025_9573_MOESM3_ESM.zip"
EXT = SUP / "MOESM3_ESM"
TAB = EXT / "2024-07-15352A-s3"
if not (TAB / "PsychAD_SupplementaryTable4.xlsx").exists():
    with zipfile.ZipFile(ZIP) as z:
        z.extractall(EXT, members=[m for m in z.namelist() if "__MACOSX" not in m])
T4, T10 = TAB / "PsychAD_SupplementaryTable4.xlsx", TAB / "PsychAD_SupplementaryTable10.csv"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison", "gene", "rank", "stat",
        "padj", "source", "selection"]
rows, unres = [], []

# PsychAD subclass -> schema broad class
BROAD = {"Astro": "Astrocytes", "Micro": "Microglia", "PVM": "Immune (other)", "Adaptive": "T cells",
         "Oligo": "Oligodendrocytes", "OPC": "OPCs", "Endo": "Endothelial", "PC": "Pericytes", "SMC": "Pericytes",
         "VLMC": "Fibroblasts"}


def broad(sub):
    if sub in BROAD:
        return BROAD[sub]
    if sub.startswith("EN_"):
        return "Excitatory neurons"
    if sub.startswith("IN_"):
        return "Inhibitory neurons"
    return "All"


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


# ---------------------------------------------------------------- DE per subclass, PD and DLBD vs control (Table 4)
t4 = pd.read_excel(T4, sheet_name=["PD", "DLBD"])
EXTRA = {"IN_LAMP5_RELN": cid(7), "IN_ADARB2": cid(7), "IN_PVALB": cid(7), "VLMC": cid(8)}
for dx, d in t4.items():
    sig = d[d.FDR < 0.05]
    for sub in sorted(d.assay.unique()):
        s = sig[sig.assay == sub]
        for sign, kind in [(1, "de_up"), (-1, "de_down")]:
            x = s[np.sign(s.estimate) == sign].sort_values("estimate", ascending=(sign < 0))
            if x.empty:
                continue
            add(f"{dx}_vs_control_{sub}_{kind}", cid(1) + EXTRA.get(sub, []), kind, broad(sub), sub,
                f"{dx} vs control", x.ID.tolist(), f"{T4.name}|{dx}",
                "dreamlet pseudobulk, meta-analysed across tissue sources: FDR<0.05 (all genes, shared and "
                "disorder-specific); ranked by |logFC| (estimate)", stat=x.estimate.tolist(), padj=x.FDR.tolist())
    for sub, cl in EXTRA.items():
        if sig[sig.assay == sub].empty:
            unresolved(cl, f"{dx}_vs_control_{sub}", f"0 {dx} DEGs at FDR<0.05 in {sub}")

# ---------------------------------------------------------------- cross-disorder shared genes (C09)
shared = pd.read_excel(T4, sheet_name="Cross-disorder Shared Genes", header=None)[0].dropna().astype(str).tolist()
pp = t4["PD"].groupby("ID").crossDis_PP.max()
s = pd.Series({g: pp.get(g, np.nan) for g in shared}).sort_values(ascending=False, kind="stable")
add("crossdisorder_shared_genes", cid(9), "pathway_genes", "All", "", "8 disorders vs control (shared component)",
    s.index.tolist(), f"{T4.name}|Cross-disorder Shared Genes",
    "all genes listed (no direction); ranked by max crossDis_PP across subclasses (PD sheet)",
    stat=s.round(4).tolist())

# ---------------------------------------------------------------- subclass / subtype markers (C03-C05)
t10 = pd.read_csv(T10, encoding="utf-8-sig")
STATES = [  # name, claims, state, cell_type, genes from the name, Table 10 filter
    ("VLMC", cid(3), "VLMC", "Fibroblasts", [], (t10.subclass == "VLMC") & t10.subtype.isna()),
    ("VLMC_DCDC2", cid(5), "VLMC_DCDC2", "Fibroblasts", ["DCDC2"], t10.subtype == "VLMC_DCDC2"),
    ("IN_LAMP5_RELN", cid(4), "IN_LAMP5_RELN", "Inhibitory neurons", ["LAMP5", "RELN"],
     (t10.subclass == "IN_LAMP5_RELN") & t10.subtype.isna()),
    ("IN_LAMP5_LHX6", cid(4), "IN_LAMP5_LHX6", "Inhibitory neurons", ["LAMP5", "LHX6"],
     (t10.subclass == "IN_LAMP5_LHX6") & t10.subtype.isna()),
    ("IN_ADARB2", cid(4), "IN_ADARB2", "Inhibitory neurons", ["ADARB2"],
     (t10.subclass == "IN_ADARB2") & t10.subtype.isna()),
]
for name, cl, st, ct, ngenes, filt in STATES:
    panel = t10.loc[filt, "marker"].tolist()
    genes = ngenes + [g for g in panel if g not in ngenes]
    add(f"{name}_markers", cl, "state_markers", ct, st, "", genes,
        "text; " + f"{T10.name}|", "text list: genes in the subclass name first, then Xenium panel markers for that "
        "subclass/subtype (Supplementary Table 10, file order); no marker statistics supplied")

# ---------------------------------------------------------------- named genes (C11)
add("MAPT_H1_upregulated_genes", cid(11), "named_genes", "All", "", "H1 vs H2 haplotype",
    ["LRRC37A2", "LRRC37A3", "PLEKHM1", "PDE8B", "HNRNPK"], "text", "genes named in the text (order as written)")

# ---------------------------------------------------------------- pathway terms named only
unresolved(cid(7), "GO BP: neuronal development; synaptic signalling (IN_LAMP5_RELN, IN_ADARB2, IN_PVALB)",
           "pathway themes named only (no term table in the supplement); per-subclass PD/DLBD DE lists provided")
unresolved(cid(8), "GO BP: vasculature development (VLMC)",
           "pathway theme named only; PD has 0 VLMC DEGs, DLBD VLMC DE lists provided")
unresolved(cid(9), "GO BP: mRNA splicing / mRNA processing; protein localisation to mitochondrion",
           "pathway themes named only; shared-gene list provided as crossdisorder_shared_genes")

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

"""Signatures for lee2023 (Lee et al. 2023, Sci Adv; SN snRNA + snATAC + bulk H3K27ac).

Reads Tables S3 (DEGs), S4 (dysregulated cREs) and S5 (ABC target genes) from the
supplement zip (extracted into <supp>/tables_s1_to_s10/) and writes
claims/signatures/lee2023.tsv, lee2023_regions.tsv and lee2023_unresolved.tsv.
Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/lee2023.py
"""
from pathlib import Path
import zipfile
import pandas as pd
import yaml

PAPER = "lee2023"
ROOT = Path(__file__).resolve().parents[2]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

SUP = cfg.RAW / "supplements" / PAPER
ZIP = SUP / "sciadv.abo2467_tables_s1_to_s10.zip"
TAB = SUP / "tables_s1_to_s10"
OUT = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison",
        "gene", "rank", "stat", "padj", "source", "selection"]
RCOLS = ["signature_id", "paper", "claim_ids", "cell_type", "comparison", "chrom", "start", "end",
         "direction", "stat", "padj", "nearest_gene", "source"]
CT = {"DopaN": "DA neurons", "GabaN": "Inhibitory neurons", "Oligo": "Oligodendrocytes", "OPC": "OPCs",
      "Ast": "Astrocytes", "Micro": "Microglia", "Endo": "Endothelial", "Peri": "Pericytes"}

if not (TAB / "abo2467_Table_S3.xlsx").exists():
    TAB.mkdir(exist_ok=True)
    zipfile.ZipFile(ZIP).extractall(TAB)

rows, rrows = [], []


def add(name, claims, kind, ct, state, comp, genes, stats, padj, source, selection):
    for i, g in enumerate(genes):
        rows.append(dict(signature_id=f"{PAPER}__{name}", paper=PAPER, claim_ids=";".join(claims), kind=kind,
                         cell_type=ct, cell_state=state, comparison=comp, gene=str(g).strip(), rank=i + 1,
                         stat="" if stats is None else stats[i], padj="" if padj is None else padj[i],
                         source=source, selection=selection))


# ---- Table S3: donor-pseudobulk edgeR DEGs per cell type ---------------------------------------
s3 = pd.read_excel(TAB / "abo2467_Table_S3.xlsx", header=1)
s3.columns = ["cell", "dir", "gene", "lfc", "p", "fdr"]
extra = {("DopaN", "Down"): ["lee2023-C03", "lee2023-C08"], ("DopaN", "Up"): ["lee2023-C06", "lee2023-C08"],
         ("Ast", "Down"): ["lee2023-C04"], ("OPC", "Down"): ["lee2023-C04"],
         ("Micro", "Down"): ["lee2023-C05", "lee2023-C09"], ("Micro", "Up"): ["lee2023-C06", "lee2023-C09"],
         ("Oligo", "Up"): ["lee2023-C07", "lee2023-C10"], ("Oligo", "Down"): ["lee2023-C10"],
         ("Endo", "Up"): ["lee2023-C07"]}
for (cell, d), g in s3.groupby(["cell", "dir"]):
    g = g.sort_values("lfc", ascending=(d == "Down"))
    kind = "de_up" if d == "Up" else "de_down"
    add(f"{cell}_PD_vs_ctrl_{d.lower()}", extra.get((cell, d), []), kind, CT[cell], "", "PD vs control",
        g.gene.tolist(), g.lfc.round(4).tolist(), g.fdr.tolist(), "abo2467_Table_S3.xlsx|Table S3",
        "authors' DEGs: edgeR LRT on donor pseudobulk, FDR<0.05 & |log2FC|>0.5; ranked by log2FC")

# ---- AGTR1 DopaN state (composition claims) and named genes --------------------------------------
add("AGTR1_DopaN_markers", ["lee2023-C01", "lee2023-C02"], "state_markers", "DA neurons", "AGTR1+ DopaN", "",
    ["AGTR1", "TH", "SLC6A3"], None, None, "text", "text list (AGTR1+ within TH+/SLC6A3+ DopaN; AGTR1- = complement)")
claims = yaml.safe_load(open(ROOT / "claims" / f"{PAPER}.yaml"))
for c in claims:
    if c["id"] in ("lee2023-C08", "lee2023-C09", "lee2023-C10", "lee2023-C18", "lee2023-C19"):
        add(f"{c['id'].split('-')[1]}_named", [c["id"]], "named_genes", c["cell_type"], "", c["comparison"],
            c["genes"], None, None, "text", "genes named in claim text")

# ---- Regions: Table S4 dysregulated cREs (bulk H3K27ac, not cell-type resolved) -------------------
s4 = pd.read_excel(TAB / "abo2467_Table_S4.xlsx", header=1)
coord = s4["cRE ID"].str.extract(r"^(chr[^:]+):(\d+)-(\d+)$")
assert coord.notna().all().all()
for lab in ("Up", "Down"):
    m = s4.Label == lab
    for (_, r), (_, cc) in zip(s4[m].iterrows(), coord[m].iterrows()):
        rrows.append(dict(signature_id=f"{PAPER}__H3K27ac_cRE_{lab.lower()}", paper=PAPER,
                          claim_ids="lee2023-C11;lee2023-C12;lee2023-C16", cell_type="All", comparison="PD vs control",
                          chrom=cc[0], start=int(cc[1]), end=int(cc[2]), direction=lab.lower(), stat=r.Log2FC,
                          padj=r.FDR, nearest_gene="",
                          source="abo2467_Table_S4.xlsx|Table S4 (bulk H3K27ac ChIP, edgeR QL FDR<0.05; stat=log2FC)"))

# ---- Regions: Table S5 ABC-linked cREs per cell type (target gene in nearest_gene) ----------------
s5 = pd.read_excel(TAB / "abo2467_Table_S5.xlsx", sheet_name=None, header=1)
rclaims = {("Oligo", "DownPeak"): "lee2023-C18", ("DopaN", "UpPeak"): "lee2023-C19",
           ("Oligo", "UpPeak"): "lee2023-C19", ("Micro", "UpPeak"): "lee2023-C19"}
for cell, d in s5.items():
    d = d[d["cRE type"].isin(["UpPeak", "DownPeak"])]
    for typ, g in d.groupby("cRE type"):
        g = g.sort_values("ABC score", ascending=False)
        cc = g["cRE ID"].str.extract(r"^(chr[^:]+):(\d+)-(\d+)$")
        for (_, r), (_, c) in zip(g.iterrows(), cc.iterrows()):
            rrows.append(dict(signature_id=f"{PAPER}__{cell}_ABC_{typ}", paper=PAPER,
                              claim_ids=rclaims.get((cell, typ), ""), cell_type=CT[cell], comparison="PD vs control",
                              chrom=c[0], start=int(c[1]), end=int(c[2]),
                              direction="up" if typ == "UpPeak" else "down", stat=r["ABC score"], padj="",
                              nearest_gene=r["Gene ID"],
                              source=f"abo2467_Table_S5.xlsx|{cell} (ABC-linked target gene, not nearest; stat=ABC score)"))

sig = pd.DataFrame(rows, columns=COLS)
reg = pd.DataFrame(rrows, columns=RCOLS)

unres = pd.DataFrame([
    ("lee2023-C03", "GO: mitochondrion organization / mitochondrial function (fig. S6B)",
     "GO term names only (fig. S6B); resolve members from MSigDB. DopaN down-DEGs attached"),
    ("lee2023-C04", "GO: regulation of neurogenesis (fig. S6B)", "GO term name only; Ast/OPC down-DEGs attached"),
    ("lee2023-C05", "GO: immune response (fig. S6B)", "GO term name only; Micro down-DEGs attached"),
    ("lee2023-C06", "GO: macroautophagy (fig. S6B)", "GO term name only; DopaN/Micro up-DEGs attached"),
    ("lee2023-C07", "GO: protein folding / protein stabilization (fig. S6B)", "GO term name only; Oligo/Endo up-DEGs attached"),
    ("lee2023-C11", "chromatin_da", "cREs are from bulk H3K27ac (no cell-type label in Table S4); regions given as All"),
    ("lee2023-C12", "chromatin_da", "uses Table S4 bulk cREs; regions given"),
    ("lee2023-C18", "GO: myelination (target genes of Oligo down cREs)", "GO term name only; Oligo DownPeak ABC regions + named genes given"),
    ("lee2023-C19", "GO: autophagy / protein folding (target genes of up cREs)", "GO term names only; ABC UpPeak regions + named genes given"),
], columns=["claim_id", "term", "reason"])

OUT.mkdir(parents=True, exist_ok=True)
sig.to_csv(OUT / f"{PAPER}.tsv", sep="\t", index=False)
reg.to_csv(OUT / f"{PAPER}_regions.tsv", sep="\t", index=False)
unres.to_csv(OUT / f"{PAPER}_unresolved.tsv", sep="\t", index=False)

assert sig.gene.str.len().gt(0).all() and (sig.groupby("signature_id")["rank"].min() == 1).all()
print(sig.groupby(["signature_id", "kind"]).size().to_string())
print(reg.groupby(["signature_id", "claim_ids"]).size().to_string())
need = {c["id"] for c in claims if c["test_template"] in ("state_signature", "geneset_shift", "de_pseudobulk")
        or (c["test_template"] == "composition_donor" and c.get("cell_state"))}
covered = set(";".join(sig.claim_ids).split(";")) | set(";".join(reg.claim_ids).split(";"))
print(f"\nrequired claims: {len(need)}; in signatures: {len(need & covered)}; MISSING:",
      sorted(need - covered - set(unres.claim_id)) or "none")
chrom = {c["id"] for c in claims if c["test_template"] == "chromatin_da"}
print("chromatin claims with regions:", sorted(chrom & covered), "neither:", sorted(chrom - covered - set(unres.claim_id)) or "none")

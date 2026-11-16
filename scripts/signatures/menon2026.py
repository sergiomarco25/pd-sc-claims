"""Signatures for menon2026 (Menon et al. 2026 preprint; multi-region multiome, Banner cohort).

media-2.xlsx: SupTable2 (DEGs), SupTable3 (Enrichr on DEGs), SupTable4 (DA peaks linked to DEGs),
SupTable5 (LB-stage trajectory DEGs). Each sheet has a text preamble; the header row follows the
blank line after 'Column Explanations'.
Writes claims/signatures/menon2026.tsv, menon2026_regions.tsv, menon2026_unresolved.tsv.
Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/menon2026.py
"""
from pathlib import Path
import pandas as pd
import yaml

PAPER = "menon2026"
ROOT = Path(__file__).resolve().parents[2]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

SUP = cfg.RAW / "supplements" / PAPER
XL = SUP / "media-2.xlsx"
OUT = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison",
        "gene", "rank", "stat", "padj", "source", "selection"]
RCOLS = ["signature_id", "paper", "claim_ids", "cell_type", "comparison", "chrom", "start", "end",
         "direction", "stat", "padj", "nearest_gene", "source"]
CT = {"Oligo": "Oligodendrocytes", "DRD1-MSNs": "MSNs", "DRD2-MSNs": "MSNs", "DRD1-DRD2-MSNs": "MSNs",
      "Excitatory": "Excitatory neurons", "RORB-Neuron": "Excitatory neurons", "Astrocyte": "Astrocytes",
      "Bergmann-Glia": "Astrocytes", "Granule": "Neurons (other)", "Dopaminergic": "DA neurons",
      "SST": "Inhibitory neurons", "LAMP5": "Inhibitory neurons", "VIP": "Inhibitory neurons",
      "PVALB": "Inhibitory neurons", "Chandelier": "Inhibitory neurons", "GABA-Put": "Inhibitory neurons",
      "GABA-SN": "Inhibitory neurons", "MLI": "Inhibitory neurons", "Golgi": "Inhibitory neurons",
      "Purkinje": "Inhibitory neurons", "Microglia": "Microglia", "OPC": "OPCs", "PV-Fibroblast": "Fibroblasts",
      "Fibroblast": "Fibroblasts", "SMC": "Pericytes", "Pericyte": "Pericytes", "Capillary": "Endothelial",
      "PV-Macrophage": "Immune (other)", "T-cell": "T cells", "Ependymal": "All"}
CEREBELLAR = {"Granule", "Purkinje", "Bergmann-Glia", "Golgi", "MLI"}
rows, rrows = [], []


def load(sheet):
    d = pd.read_excel(XL, sheet_name=sheet, header=None)
    c0 = d.iloc[:, 0].tolist()
    ce = [i for i, v in enumerate(c0) if str(v).strip() == "Column Explanations"][0]
    h = [i for i in range(ce, len(c0)) if pd.isna(c0[i])][0] + 1
    t = d.iloc[h + 1:].copy()
    t.columns = [str(x) for x in d.iloc[h]]
    return t.dropna(how="all").reset_index(drop=True)


def add(name, claims, kind, ct, state, comp, genes, stats, padj, source, selection):
    if stats is None and padj is None or kind == "pathway_genes":
        genes = list(dict.fromkeys(str(g).strip() for g in genes))  # de-duplicate text/pathway lists, keep order
    for i, g in enumerate(genes):
        rows.append(dict(signature_id=f"{PAPER}__{name}", paper=PAPER, claim_ids=";".join(claims), kind=kind,
                         cell_type=ct, cell_state=state, comparison=comp, gene=str(g).strip(), rank=i + 1,
                         stat="" if stats is None else stats[i], padj="" if padj is None else padj[i],
                         source=source, selection=selection))


# ---- SupTable2: pseudobulk edgeR DEGs, PD vs control (regions pooled; no region column) -----------------
t2 = load("SupTable2 DE_GEX")
for c in ("avg_logFC", "p_val_adj"):
    t2[c] = pd.to_numeric(t2[c])
DE_CLAIMS = {"Oligo": ["menon2026-C02", "menon2026-C08", "menon2026-C09"],
             "DRD1-MSNs": ["menon2026-C03"], "DRD2-MSNs": ["menon2026-C03"], "DRD1-DRD2-MSNs": ["menon2026-C03"],
             "Excitatory": ["menon2026-C04", "menon2026-C09"],
             "Astrocyte": ["menon2026-C05", "menon2026-C08", "menon2026-C09"],
             "Dopaminergic": ["menon2026-C06"]}
for ct, g in t2.groupby("cell_type"):
    cl = DE_CLAIMS.get(ct, []) + (["menon2026-C10"] if ct in CEREBELLAR else [])
    for d, sub in (("up", g[g.avg_logFC > 0].sort_values("avg_logFC", ascending=False)),
                   ("down", g[g.avg_logFC < 0].sort_values("avg_logFC"))):
        if sub.empty:
            continue
        add(f"{ct}_PD_vs_ctrl_{d}", cl, f"de_{d}", CT[ct], ct, "PD vs control", sub.gene.tolist(),
            sub.avg_logFC.round(4).tolist(), sub.p_val_adj.tolist(), f"{XL.name}|SupTable2 DE_GEX",
            "authors' DEGs: edgeR glmLRT on pseudobulk, p_val_adj<0.05 (as listed), split by sign, ranked by logFC")

# ---- SupTable3: DA-neuron dopamine/catecholamine terms (member genes given) --------------------------------
t3 = load("SupTable3 PathwayEnrichmentDEGs")
da = t3[(t3.cell_type == "Dopaminergic") & t3.Term.str.contains("Dopamine Transport|Dopamine Uptake \\(|Catecholamine Transport|Catecholamine Uptake", regex=True)].drop_duplicates("Term")
for _, r in da.iterrows():
    gid = r.Term.split("(")[-1].rstrip(")").replace(":", "")
    genes = [x for x in str(r.Genes).split(";") if x]
    add(f"DA_{gid}", ["menon2026-C06"], "pathway_genes", "DA neurons", "", "PD vs control", genes, None,
        [r["Adjusted.P.value"]] * len(genes), f"{XL.name}|SupTable3 PathwayEnrichmentDEGs",
        f"Enrichr {r.Term} (adj.P={float(r['Adjusted.P.value']):.3g}); overlap genes with DA DEGs as listed")

# ---- SupTable5: LB-stage trajectory genes (2 unlabelled test rows per gene; take min padj) -------------------
t5 = load("SupTable5 TrajectoryDEGs")
for c in ("avg_logFC", "p_val_adj"):
    t5[c] = pd.to_numeric(t5[c])
best = t5.sort_values("p_val_adj").groupby(["cell_type", "gene"], as_index=False).first()
best = best[best.p_val_adj < 0.05]
TRAJ = {("monotonic_increase", "up"): "de_up", ("monotonic_decrease", "down"): "de_down",
        ("late_response", "up"): "de_up", ("late_response", "down"): "de_down",
        ("early_response", "up"): "de_up", ("early_response", "down"): "de_down"}
for (ct, pat, d), g in best.groupby(["cell_type", "pattern", "direction"]):
    if (pat, d) not in TRAJ:
        continue
    cl = []
    if "MSN" in ct and pat == "monotonic_increase":
        cl = ["menon2026-C17"]
    if ct == "Astrocyte" and pat == "late_response" and d == "up":
        cl = ["menon2026-C15"]
    g = g.sort_values("p_val_adj")
    add(f"{ct}_LBstage_{pat}_{d}", cl, TRAJ[(pat, d)], CT[ct], ct, f"across Unified LB stage ({pat})",
        g.gene.tolist(), g.avg_logFC.round(4).tolist(), g.p_val_adj.tolist(), f"{XL.name}|SupTable5 TrajectoryDEGs",
        "min p_val_adj<0.05 over the two test rows per gene; authors' pattern/direction; ranked by padj")

# ---- named genes / text-list states ------------------------------------------------------------------------
add("AGTR1_DA_subcluster6_markers", ["menon2026-C01"], "state_markers", "DA neurons", "AGTR1+ DA subcluster 6", "",
    ["AGTR1"], None, None, "text", "text list (Fig. 1g; subcluster marker table not in supplement)")
add("eccentric_MSN_markers", ["menon2026-C21"], "state_markers", "MSNs", "eccentric DRD1/DRD2 MSNs", "",
    ["DRD1", "DRD2"], None, None, "text", "text list (co-expression of DRD1 and DRD2)")
claims = yaml.safe_load(open(ROOT / "claims" / f"{PAPER}.yaml"))
for c in claims:
    if c["id"] in ("menon2026-C14", "menon2026-C15", "menon2026-C18"):
        add(f"{c['id'].split('-')[1]}_named", [c["id"]], "named_genes", c["cell_type"], "", c["comparison"],
            c["genes"], None, None, "text", "genes named in claim text")

# ---- regions: SupTable4 DA peaks linked to DEGs -----------------------------------------------------------
t4 = load("SupTable 4 DEG-linked Peak DA")
cc = t4.peak_coord.str.extract(r"^(chr[^:]+):(\d+)-(\d+)$")
assert cc.notna().all().all()
t4 = pd.concat([t4, cc.set_axis(["chrom", "start", "end"], axis=1)], axis=1)
t4["log2FoldChange"] = pd.to_numeric(t4.log2FoldChange)
for (cell, up), g in t4.groupby(["cell", t4.log2FoldChange > 0]):
    d = "up" if up else "down"
    g = g.sort_values("log2FoldChange", ascending=not up)
    for _, r in g.iterrows():
        rrows.append(dict(signature_id=f"{PAPER}__{cell}_DEGlinked_DApeaks_{d}", paper=PAPER, claim_ids="menon2026-C07",
                          cell_type=CT[cell], comparison="PD vs control", chrom=r.chrom, start=int(r.start), end=int(r.end),
                          direction=d, stat=r.log2FoldChange, padj=r.padj_gene, nearest_gene=r.target_gene,
                          source=f"{XL.name}|SupTable 4 (DESeq2; padj = gene-level padj_gene; nearest_gene = authors' target gene)"))

sig = pd.DataFrame(rows, columns=COLS)
reg = pd.DataFrame(rrows, columns=RCOLS)
unres = pd.DataFrame([
    ("menon2026-C08", "region-resolved DEGs (SN vs other regions)", "SupTable2 has no region column; only region-pooled DEGs given (attached)"),
    ("menon2026-C09", "region-resolved DEG overlap", "SupTable2 has no region column; region-pooled DEGs attached"),
    ("menon2026-C10", "cerebellum DEGs", "no region column; cerebellum-specific cell types (Granule, Purkinje, Bergmann-Glia, Golgi, MLI) attached as proxy"),
    ("menon2026-C16", "chromatin_da (SERPINH1 astrocyte peaks across LB stage)", "no SERPINH1 peak in SupTable4 and no stage-wise DA peak table in supplement"),
    ("menon2026-C01", "AGTR1+ DA subcluster 6 marker table", "subcluster markers not in supplement; AGTR1 text list used"),
    ("menon2026-C21", "eccentric MSN markers (Supp. Fig. 4d-f)", "marker table not in supplement; DRD1/DRD2 text list used"),
], columns=["claim_id", "term", "reason"])

OUT.mkdir(parents=True, exist_ok=True)
sig.to_csv(OUT / f"{PAPER}.tsv", sep="\t", index=False)
reg.to_csv(OUT / f"{PAPER}_regions.tsv", sep="\t", index=False)
unres.to_csv(OUT / f"{PAPER}_unresolved.tsv", sep="\t", index=False)

assert sig.gene.str.len().gt(0).all() and (sig.groupby("signature_id")["rank"].min() == 1).all()
print(sig.groupby(["signature_id", "kind"]).size().to_string())
print(reg.groupby("signature_id").size().to_string())
need = {c["id"] for c in claims if c["test_template"] in ("state_signature", "geneset_shift", "de_pseudobulk")
        or (c["test_template"] == "composition_donor" and c.get("cell_state"))}
covered = set(";".join(sig.claim_ids).split(";")) | set(";".join(reg.claim_ids).split(";"))
print(f"\nrequired claims: {len(need)}; in signatures: {len(need & covered)}; MISSING:",
      sorted(need - covered - set(unres.claim_id)) or "none")
chrom = {c["id"] for c in claims if c["test_template"] == "chromatin_da"}
print("chromatin claims with regions:", sorted(chrom & covered), "unresolved:", sorted(chrom & set(unres.claim_id)),
      "neither:", sorted(chrom - covered - set(unres.claim_id)) or "none")

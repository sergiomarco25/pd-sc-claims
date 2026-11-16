"""Signatures for clarence2025 (Clarence et al. 2025, medRxiv; AMP-PD multi-region snRNA).

The only supplementary file available is NIHPP2025.07.30.25332436V1-supplement-1.pdf
(13 pages of supplementary figures, no gene tables). Supplementary Tables 4-12 (dreamlet DEGs,
GO results, GWAS-gene DE, TF activities, LR pairs, Hotspot modules, subtype composition) are
not in the supplement folder, so only genes named in the claim text are emitted
(read from claims/clarence2025.yaml `genes`, source = text). Everything else goes to unresolved.
Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/clarence2025.py
"""
from pathlib import Path
import pandas as pd
import yaml

PAPER = "clarence2025"
ROOT = Path(__file__).resolve().parents[2]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

SUP = cfg.RAW / "supplements" / PAPER
OUT = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison",
        "gene", "rank", "stat", "padj", "source", "selection"]

xlsx = sorted(SUP.rglob("*.xlsx")) + sorted(SUP.rglob("*.csv"))
print("tabular supplementary files found:", [p.name for p in xlsx] or "none")

claims = {c["id"]: c for c in yaml.safe_load(open(ROOT / "claims" / f"{PAPER}.yaml"))}
rows = []


def add(name, cids, kind, state, genes, selection):
    c = claims[cids[0]]
    for i, g in enumerate(dict.fromkeys(genes)):
        rows.append(dict(signature_id=f"{PAPER}__{name}", paper=PAPER, claim_ids=";".join(cids), kind=kind,
                         cell_type=c["cell_type"], cell_state=state, comparison=c["comparison"] if kind != "state_markers" else "",
                         gene=g, rank=i + 1, stat="", padj="", source="text", selection=selection))


# named genes (DE claims)
add("C02_HSP_named", ["clarence2025-C02"], "named_genes", "", claims["clarence2025-C02"]["genes"], "genes named in claim text")
add("C12_GWAS_DE_named", ["clarence2025-C12"], "named_genes", "", claims["clarence2025-C12"]["genes"],
    "genes named in claim text (subset of Supp. Table 6 84 GWAS genes)")
# state text lists for composition claims
add("homeostatic_microglia_markers", ["clarence2025-C20"], "state_markers", "homeostatic microglia", ["P2RY12", "CX3CR1"], "text list")
add("proliferative_microglia_markers", ["clarence2025-C21", "clarence2025-C25"], "state_markers",
    "adaptive and proliferative microglia (Prolif_MKI67)", ["MKI67"], "text list (subtype name gene)")
add("Adapt_TMEM163_AIF1_markers", ["clarence2025-C23"], "state_markers", "Adapt_TMEM163, Adapt_AIF1", ["TMEM163", "AIF1"], "text list (subtype name genes)")
add("Adapt_IFI44L_CCL3_markers", ["clarence2025-C24"], "state_markers", "Adapt_IFI44L, Adapt_CCL3", ["IFI44L", "CCL3"], "text list (subtype name genes)")
add("Adapt_HSPA1A_markers", ["clarence2025-C25"], "state_markers", "Adapt_HSPA1A", ["HSPA1A"], "text list (subtype name gene)")
add("PVM_CD163_markers", ["clarence2025-C26"], "state_markers", "PVM_CD163", ["CD163"], "text list (subtype name gene; no DAM genes named)")

sig = pd.DataFrame(rows, columns=COLS)
NOTAB = "Supplementary Tables not in supplement folder (only a figures PDF)"
unres = pd.DataFrame([
    ("clarence2025-C01", "dreamlet DEGs per cell type (Supp. Table 4)", NOTAB),
    ("clarence2025-C02", "full HSP DE lists for mural/endothelial/OPC (Supp. Table 4)", NOTAB + "; named genes emitted"),
    ("clarence2025-C03", "GO:0006986 response to unfolded protein", "pathway name only; resolve from MSigDB/GO"),
    ("clarence2025-C04", "GO:0006974 DNA damage response", "pathway name only; resolve from MSigDB/GO"),
    ("clarence2025-C05", "GO:0090083 regulation of inclusion body assembly; autophagy", "pathway name only; resolve from MSigDB/GO"),
    ("clarence2025-C06", "CollecTRI regulon: HSF1", "TF activity; resolve regulon centrally (decoupler CollecTRI)"),
    ("clarence2025-C07", "CollecTRI regulon: FOXO3", "TF activity; resolve regulon centrally (decoupler CollecTRI)"),
    ("clarence2025-C08", "CollecTRI regulon: MYC", "TF activity; resolve regulon centrally (decoupler CollecTRI)"),
    ("clarence2025-C09", "CollecTRI regulon: AR", "TF activity; resolve regulon centrally (decoupler CollecTRI)"),
    ("clarence2025-C12", "84 GWAS genes with DE t-statistics (Supp. Table 6)", NOTAB + "; named genes emitted"),
    ("clarence2025-C19", "EN_L6_CT, EN_L6B subtype markers", NOTAB + "; no marker genes named in text"),
    ("clarence2025-C20", "microglia subtype markers (Supp. Table 12 / Fig. 7b)", NOTAB + "; text list emitted"),
    ("clarence2025-C21", "microglia subtype markers", NOTAB + "; text list emitted"),
    ("clarence2025-C23", "microglia subtype markers", NOTAB + "; subtype-name genes emitted"),
    ("clarence2025-C24", "microglia subtype markers", NOTAB + "; subtype-name genes emitted"),
    ("clarence2025-C25", "microglia subtype markers", NOTAB + "; subtype-name genes emitted"),
    ("clarence2025-C26", "DAM markers and PVM_CD163 markers", NOTAB + "; only CD163 emitted, DAM unresolved"),
], columns=["claim_id", "term", "reason"])

OUT.mkdir(parents=True, exist_ok=True)
sig.to_csv(OUT / f"{PAPER}.tsv", sep="\t", index=False)
unres.to_csv(OUT / f"{PAPER}_unresolved.tsv", sep="\t", index=False)

assert sig.gene.str.len().gt(0).all() and (sig.groupby("signature_id")["rank"].min() == 1).all()
print(sig.groupby(["signature_id", "kind"]).size().to_string())
need = {c["id"] for c in claims.values() if c["test_template"] in ("state_signature", "geneset_shift", "de_pseudobulk")
        or (c["test_template"] == "composition_donor" and c.get("cell_state"))}
covered = set(";".join(sig.claim_ids).split(";"))
print(f"\nrequired claims: {len(need)}; in signatures: {len(need & covered)}; unresolved: {len(need & set(unres.claim_id))}; MISSING:",
      sorted(need - covered - set(unres.claim_id)) or "none")
chrom = {c["id"] for c in claims.values() if c["test_template"] == "chromatin_da"}
print("chromatin claims:", sorted(chrom) or "none")

"""Signatures for zhu2024 (Zhu et al. 2024, Sci Transl Med; PFC BA9 snRNA, 6 PD / 6 control).

The supplementary tables (Science, NIHMS2101876) could not be obtained. This script therefore emits
named_genes signatures from the claim text and lists everything else as unresolved with reason
"supplement not obtainable".

If an xlsx appears later in the supplement folder (expected:
NIHMS2101876-supplement-Supplemental_Tables.xlsx; any *.xlsx is tried), the script automatically
tries to parse DEG tables from it (table S3: cell-type DEGs, PD vs control): it looks in every sheet
whose name or first rows mention S3 / DEG / differential for a header row with a gene column and a
logFC column (plus optional adjusted-p and cell-type/cluster columns), and emits de_up / de_down
signatures per cell type. Claims covered this way are dropped from the unresolved file. The parsed
output is marked "auto-parsed" in `selection`: check it by eye the first time it runs.

Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/zhu2024.py [supplement_dir]
(the optional argument overrides the supplement folder, for testing only)
"""
from pathlib import Path
import re
import sys
import pandas as pd
import yaml

PAPER = "zhu2024"
ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pdrep import config as cfg  # noqa: E402
SUP = Path(sys.argv[1]) if len(sys.argv) > 1 else \
    cfg.RAW / "supplements" / PAPER
EXPECTED = "NIHMS2101876-supplement-Supplemental_Tables.xlsx"
OUT = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison",
        "gene", "rank", "stat", "padj", "source", "selection"]
C = lambda *n: [f"{PAPER}-C{int(i):02d}" for i in n]  # noqa: E731
COMP = "PD vs control"
NO_SUP = "supplement not obtainable"

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


# ---------------- named genes from the claim text ----------------
NAMED = [
    ("Tcell_markers_named", C(1), "state_markers", "T cells", "brain-resident T cells", "", ["CD3E", "SKAP1", "GNLY", "CD96"]),
    ("Neuron_axonogenesis_up_named", C(5), "named_genes", "Neurons (other)", "ExN, InN", COMP,
     ["EPHA6", "ROBO1", "GAB1", "PARD3", "NLGN4Y", "CBLN2"]),
    ("Astro_metal_up_named", C(6), "named_genes", "Astrocytes", "", COMP, ["MT1G", "MT1F", "MT3"]),
    ("Astro_coagulation_down_named", C(7), "named_genes", "Astrocytes", "", COMP, ["VAV3", "LYN", "UBASH3B"]),
    ("Micro_activation_up_named", C(8), "named_genes", "Microglia", "", COMP,
     ["FOS", "STAT3", "RUNX1", "BCL6", "IL4R", "JAK3", "NFKBIZ", "IL15", "NFKBID"]),
    ("Oligo_up_named", C(9), "named_genes", "Oligodendrocytes", "", COMP,
     ["HSPA1A", "HSPA1B", "FMN1", "P2RX7", "RAPGEF3", "FCHSD2"]),
    ("OPC_up_named", C(10), "named_genes", "OPCs", "", COMP, ["VEGFA", "HSPA1A", "HSPA1B", "ZNF365", "SYT1", "MAPT"]),
    ("Tcell_up_named", C(11), "named_genes", "T cells", "", COMP, ["IL7R", "IL15", "MALT1", "CAMK4", "CD28", "ZAP70"]),
    ("ExN_chaperone_down_named", C(12), "named_genes", "Excitatory neurons", "", COMP,
     ["HSPB1", "HSP90AB1", "DNAJB6", "HSPA1A", "HSPH1", "DNAJA4"]),
    ("InN_UPR_down_named", C(13), "named_genes", "Inhibitory neurons", "", COMP, ["HSPD1", "HSP90AA1", "HSPH1", "PLCG2"]),
    ("PD_genes_up_named", C(16), "named_genes", "All", "", COMP, ["SNCA", "LRRK2"]),
    # authors write GBA (current HGNC: GBA1) and 'PRNK' (typo for PRKN)
    ("ExN_PD_genes_down_named", C(17), "named_genes", "Excitatory neurons", "", COMP, ["PRKN", "PINK1", "GBA", "PARK7"]),
    ("Oligo_UTMOST_up_named", C(19), "named_genes", "Oligodendrocytes", "", COMP, ["MAPT", "TMEM163", "KAT8"]),
    ("ExN_UTMOST_down_named", C(20), "named_genes", "Excitatory neurons", "", COMP,
     ["LRRC37A2", "BCL7C", "RNF40", "DCTPP1", "LRRC37A"]),
    ("Astro_shared_AD_up_named", C(25), "named_genes", "Astrocytes", "shared PD/AD DEGs", COMP,
     ["MT2A", "MT1E", "MRAS", "COL27A1", "PLEC", "MT1M", "MT1G"]),
    ("Astro_shared_AD_down_named", C(25), "named_genes", "Astrocytes", "shared PD/AD DEGs", COMP, ["ALDH1A1", "BAALC"]),
    ("Micro_shared_AD_up_named", C(25), "named_genes", "Microglia", "shared PD/AD DEGs", COMP,
     ["SPP1", "APOE", "APOC1", "CACNA1A", "IRF8"]),
    ("Micro_shared_AD_down_named", C(25), "named_genes", "Microglia", "shared PD/AD DEGs", COMP, ["SYTL3", "IFI44L"]),
    ("Other_shared_AD_named", C(25), "named_genes", "All", "shared PD/AD DEGs (cell type / direction not resolved)", COMP,
     ["SLC38A2", "CA2", "UTY", "RASGRF2"]),
]
for name, cl, kind, ct, st, comp, genes in NAMED:
    add(name, cl, kind, ct, st, comp, genes, None, None, "text", "text list")

# ---------------- optional: auto-parse DEG tables if the supplement appears ----------------
CT_PATTERNS = [(r"^ex|excit", "Excitatory neurons"), (r"^in(?!t)|inhib", "Inhibitory neurons"),
               (r"astro", "Astrocytes"), (r"^mg|micro", "Microglia"), (r"opc|precursor", "OPCs"),
               (r"oligo|^olg|^odc", "Oligodendrocytes"), (r"endo", "Endothelial"), (r"^t[ _-]?cell|^t$", "T cells")]
CT_CLAIMS = {"Excitatory neurons": C(4, 5, 12, 16, 17, 20), "Inhibitory neurons": C(4, 5, 13, 16),
             "Astrocytes": C(6, 7, 16, 25), "Microglia": C(8, 16, 25), "Oligodendrocytes": C(9, 16, 19, 25),
             "OPCs": C(10, 16), "T cells": C(11, 16), "Endothelial": C(16)}


def map_ct(label):
    s = str(label).strip().lower()
    for pat, ct in CT_PATTERNS:
        if re.search(pat, s):
            return ct
    return None


def find_col(cols, pats):
    for c in cols:
        if any(re.search(p, c, re.I) for p in pats):
            return c
    return None


def parse_de(xlsx):
    covered = set()
    xl = pd.ExcelFile(xlsx)
    for sh in xl.sheet_names:
        raw = xl.parse(sh, header=None, nrows=None)
        head_txt = " ".join(raw.head(5).astype(str).values.ravel())
        if not re.search(r"S3\b|DEG|differential", sh + " " + head_txt, re.I):
            continue
        hdr = None
        for i in range(min(20, len(raw))):
            vals = [str(v).strip() for v in raw.iloc[i].values]
            if find_col(vals, [r"^gene", r"symbol", r"^features?$"]) and find_col(vals, [r"log.?2?.?fc", r"logfold"]):
                hdr = i
                break
        if hdr is None:
            print(f"  [auto-parse] sheet '{sh}': no gene/logFC header found, skipped")
            continue
        d = raw.iloc[hdr + 1:].copy()
        d.columns = [str(v).strip() for v in raw.iloc[hdr].values]
        gcol = find_col(d.columns, [r"^gene", r"symbol", r"^features?$"])
        fcol = find_col(d.columns, [r"log.?2?.?fc", r"logfold"])
        pcol = find_col(d.columns, [r"adj", r"fdr", r"q.?val"])
        ccol = find_col(d.columns, [r"^cluster", r"cell.?type", r"^celltype", r"^cell$", r"^ident"])
        d = d[d[gcol].notna()]
        d[fcol] = pd.to_numeric(d[fcol], errors="coerce")
        if pcol:
            d[pcol] = pd.to_numeric(d[pcol], errors="coerce")
        groups = d.groupby(ccol) if ccol else [(sh, d)]
        for label, g in groups:
            ct = map_ct(label)
            if ct is None:
                print(f"  [auto-parse] sheet '{sh}': cell type '{label}' not mapped, skipped")
                continue
            if pcol and ct != "T cells":  # authors: T-cell DEGs were selected by |logFC|>0.75 without p filter
                g = g[g[pcol] < 0.05]
            g = g[g[fcol].abs() > (0.75 if ct == "T cells" else 0.25)]
            lab = re.sub(r"[^A-Za-z0-9]+", "_", str(label)).strip("_")
            for kind, sub in (("de_up", g[g[fcol] > 0]), ("de_down", g[g[fcol] < 0])):
                sub = sub.reindex(sub[fcol].abs().sort_values(ascending=False).index).drop_duplicates(gcol)
                if len(sub):
                    add(f"{lab}_PD_vs_control_{kind.split('_')[1]}", CT_CLAIMS[ct], kind, ct, str(label).strip(), COMP,
                        sub[gcol], sub[fcol].round(4).tolist(),
                        sub[pcol].tolist() if pcol else None, f"{Path(xlsx).name}|{sh}",
                        f"auto-parsed: {'padj<0.05 & ' if pcol and ct != 'T cells' else ''}|{fcol}|>"
                        f"{0.75 if ct == 'T cells' else 0.25}; ranked by |{fcol}| (verify parsing)")
                    covered.update(CT_CLAIMS[ct])
    return covered


xlsxs = sorted(SUP.glob("*.xlsx"), key=lambda p: (p.name != EXPECTED, p.name)) if SUP.exists() else []
covered_by_table = set()
for x in xlsxs:
    print(f"supplement found: {x.name}; trying to parse DEG tables")
    covered_by_table |= parse_de(x)
if not xlsxs:
    print(f"no supplement xlsx in {SUP}: emitting text-named genes only")

# ---------------- unresolved ----------------
TABLE_CLAIMS = {
    1: "brain-resident T cell markers beyond the 4 named genes (fig. S1/S3)",
    4: "table S3 cell-type DEGs (ExN / InN up and down lists)",
    5: "table S3 neuron DEGs", 6: "table S3 astrocyte DEGs", 7: "table S3 astrocyte DEGs",
    8: "table S3 microglia DEGs", 9: "table S3 oligodendrocyte DEGs", 10: "table S3 OPC DEGs",
    11: "table S3 T cell DEGs", 12: "table S3 excitatory neuron DEGs", 13: "table S3 inhibitory neuron DEGs",
    16: "table S3 DEGs (cell types with SNCA up not named in text)", 17: "table S3 excitatory neuron DEGs",
    19: "UTMOST PD GWAS gene list (Fig. 3C) and oligodendrocyte DEGs", 20: "UTMOST PD GWAS gene list (Fig. 3E) and ExN DEGs",
    25: "PD/AD shared DEG list (Fig. 6A)",
}
for n, term in TABLE_CLAIMS.items():
    cid = C(n)[0]
    if cid in covered_by_table and n not in (1, 19, 20, 25):
        continue
    unres(cid, term, NO_SUP + ("; named genes from text emitted" if n != 4 else "; no genes named in text"))
unres(C(12)[0], "GO BP: protein folding / response to unfolded protein (Fig. 1H, 2C)",
      "pathway named only; resolve from MSigDB (GO:0006457 protein folding, GO:0006986 response to unfolded protein)")

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

claims = yaml.safe_load(open(ROOT / "claims" / f"{PAPER}.yaml"))
need = {c["id"] for c in claims if c["test_template"] in ("state_signature", "geneset_shift", "de_pseudobulk")
        or (c["test_template"].startswith("composition") and (c.get("cell_state") or "").strip())}
in_sig = {x for s in sig.claim_ids.unique() for x in s.split(";")}
in_un = set(un.claim_id)
print(f"\ncoverage: {len(need)} claims need signatures; {len(need & in_sig)} in signatures; "
      f"{len(need & in_un - in_sig)} unresolved only; missing: {sorted(need - in_sig - in_un)}")

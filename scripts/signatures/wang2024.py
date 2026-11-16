"""Gene signatures for wang2024 claims (see claims/SIGNATURES.md).

Source: sciadv.adi8287_tables_s1_to_s9.zip -> "adi8287_Tables S1 to S9.xlsx" (extracted into
raw/supplements/wang2024/extracted/). Sheets used:
  S2.Cluster_markers        cluster markers (Seurat, FC>=1.2, Bonferroni padj<=0.05)
  S5.PD-vs-Control DEGs     per-cluster PD-vs-control DEGs (MAST, FC>=1.2, Bonferroni padj<=0.05)
  S6. MSigDB of PD DEGs     pathway enrichment of the DEGs, with overlap genes per term

Cluster -> broad class was assigned from the S2 top markers (and the paper's naming):
C0/C3 oligodendrocytes, C1 microglia, C2 astrocytes, C4 endothelial, C5 OPCs, C6_* neurons,
C6_2 and C7_3 DA neurons, C7_1 GABAergic, C7_2 glutamatergic, C8 pericytes, C9 RIT2-enriched
neurons, C10 fibroblasts, C11 T cells.

Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/wang2024.py
"""
import zipfile
from pathlib import Path

import pandas as pd
import yaml

PAPER = "wang2024"
ROOT = Path(__file__).resolve().parents[2]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

SUPP = cfg.RAW / "supplements" / PAPER
ZIP = SUPP / "sciadv.adi8287_tables_s1_to_s9.zip"
XLSX = SUPP / "extracted" / "adi8287_Tables S1 to S9.xlsx"
OUT = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison",
        "gene", "rank", "stat", "padj", "source", "selection"]

if not XLSX.exists():
    XLSX.parent.mkdir(exist_ok=True)
    zipfile.ZipFile(ZIP).extract(XLSX.name, XLSX.parent)

CELL_TYPE = {"C0": "Oligodendrocytes", "C1": "Microglia", "C2": "Astrocytes", "C3": "Oligodendrocytes",
             "C4": "Endothelial", "C5": "OPCs", "C6": "Neurons (other)", "C6_0": "Neurons (other)",
             "C6_1": "Neurons (other)", "C6_2": "DA neurons", "C6_3": "Neurons (other)",
             "C6_4": "Neurons (other)", "C6_5": "Neurons (other)", "C7_0": "Neurons (other)",
             "C7_1": "Inhibitory neurons", "C7_2": "Excitatory neurons", "C7_3": "DA neurons",
             "C8": "Pericytes", "C9": "Neurons (other)", "C10": "Fibroblasts", "C11": "T cells"}
STATE = {"C9": "c9 RIT2-enriched neurons", "C6_2": "c6_2", "C7_3": "c7_3"}
NEURONAL = lambda c: c.startswith("C6") or c.startswith("C7") or c == "C9"

# ---------------------------------------------------------------- helpers
rows, unresolved = [], []


def add_table(sig, claim_ids, kind, cell_type, cell_state, comparison, df, gene_col, stat_col,
              padj_col, source, selection):
    df = df.reset_index(drop=True)
    for i, r in df.iterrows():
        rows.append(dict(signature_id=f"{PAPER}__{sig}", paper=PAPER, claim_ids=";".join(claim_ids),
                         kind=kind, cell_type=cell_type, cell_state=cell_state, comparison=comparison,
                         gene=str(r[gene_col]).strip(), rank=i + 1,
                         stat=r[stat_col] if stat_col else "", padj=r[padj_col] if padj_col else "",
                         source=source, selection=selection))


def add_text(sig, claim_ids, kind, cell_type, cell_state, comparison, genes, selection="text list"):
    add_table(sig, claim_ids, kind, cell_type, cell_state, comparison,
              pd.DataFrame({"gene": genes}), "gene", None, None, "text", selection)


def unres(claim_id, term, reason):
    unresolved.append(dict(claim_id=claim_id, term=term, reason=reason))


def sheet(name):
    d = pd.read_excel(XLSX, sheet_name=name, header=2)
    return d.dropna(how="all")


# ---------------------------------------------------------------- state markers (S2)
mk = sheet("S2.Cluster_markers")
SRC = f"{XLSX.name}|S2.Cluster_markers"
marker_claims = {"C9": ["wang2024-C01", "wang2024-C02", "wang2024-C03", "wang2024-C23"],
                 "C6_2": ["wang2024-C04"], "C7_3": ["wang2024-C05", "wang2024-C18"]}
for cl, cids in marker_claims.items():
    d = mk[(mk.Cluster == cl) & (mk.p_val_adj < 0.05) & (mk.avg_log2FC > 0)]
    d = d.sort_values("avg_log2FC", ascending=False).head(100)
    add_table(f"{cl}_markers", cids, "state_markers", CELL_TYPE[cl], STATE[cl], "", d, "Symbol",
              "avg_log2FC", "p_val_adj", SRC,
              "p_val_adj<0.05 (Bonferroni) & avg_log2FC>0 vs all other clusters, top 100 by avg_log2FC")

# ---------------------------------------------------------------- PD-vs-control DEGs (S5)
de = sheet("S5.PD-vs-Control DEGs")
SRC = f"{XLSX.name}|S5.PD-vs-Control DEGs"
ALL = ["wang2024-C08", "wang2024-C09", "wang2024-C10", "wang2024-C15"]
per_cluster = {"C9": ["wang2024-C07", "wang2024-C13", "wang2024-C20", "wang2024-C21", "wang2024-C25"],
               "C7_3": ["wang2024-C12", "wang2024-C13", "wang2024-C20", "wang2024-C24"],
               "C6_2": ["wang2024-C13"], "C1": ["wang2024-C14", "wang2024-C22"],
               "C0": ["wang2024-C22"], "C3": ["wang2024-C22"]}
for cl in sorted(de.Cluster.unique(), key=lambda c: [int(x) for x in c[1:].split("_")]):
    cids = ALL + (["wang2024-C11", "wang2024-C24"] if NEURONAL(cl) else []) + per_cluster.get(cl, [])
    cids = list(dict.fromkeys(cids))
    for kind, sign in (("de_up", 1), ("de_down", -1)):
        d = de[(de.Cluster == cl) & (de.p_val_adj <= 0.05) & (sign * de.avg_log2FC > 0)]
        d = d.sort_values("avg_log2FC", ascending=(sign < 0))
        if d.empty:
            continue
        add_table(f"{cl}_PDvsCtrl_{kind[3:]}", cids, kind, CELL_TYPE[cl], STATE.get(cl, cl.lower()),
                  "PD vs control", d, "Symbol", "avg_log2FC", "p_val_adj", SRC,
                  "all DEGs in Table S5 (MAST, Bonferroni p_val_adj<=0.05, FC>=1.2), ranked by avg_log2FC")

# ---------------------------------------------------------------- pathway genes (S6 overlap genes)
pw = sheet("S6. MSigDB of PD DEGs")
SRC = f"{XLSX.name}|S6. MSigDB of PD DEGs"


def add_pathway(sig, cids, cell_type, cell_state, terms_regex, direction, clusters, desc):
    d = pw[pw["MSigDB Terms"].str.contains(terms_regex) & (pw["Direction of change"] == direction)
           & pw.Cluster.map(clusters)]
    genes = d.Genes.str.split(";").explode().str.strip()
    cnt = genes.value_counts().rename_axis("gene").reset_index(name="n")
    cnt = cnt.sort_values(["n", "gene"], ascending=[False, True])
    terms = ",".join(sorted(d["MSigDB Terms"].unique()))
    add_table(sig, cids, "pathway_genes", cell_type, cell_state, "PD vs control", cnt, "gene", "n", None, SRC,
              f"{desc}; union of overlap genes of {direction} terms [{terms}] in clusters "
              f"[{','.join(sorted(d.Cluster.unique()))}]; stat = n cluster-term hits")


add_pathway("ribosome_translation_up", ["wang2024-C08", "wang2024-C12"], "All", "",
            r"RIBOSOM|TRANSLATION_ELONGATION|TRANSLATION_INITIATION|^REACTOME_TRANSLATION$|SRP_DEPENDENT", "Up",
            lambda c: True, "ribosome/translation terms")
add_pathway("metallothionein_up", ["wang2024-C09"], "All", "", r"METALLOTHIONEIN|RESPONSE_TO_METAL_IONS|ZINC_HOMEOSTASIS",
            "Up", lambda c: True, "metallothionein / metal-ion terms")
add_pathway("heat_shock_up", ["wang2024-C10"], "All", "", r"HEAT_STRESS|HSF1", "Up", lambda c: True,
            "heat-shock / HSF1 terms")
add_pathway("neuronal_synaptic_vesicle_down", ["wang2024-C11"], "Neurons (other)", "", r"SYNAP|VESICLE", "Down",
            NEURONAL, "synapse / vesicle terms, neuronal clusters only")
add_pathway("C9_synaptic_down", ["wang2024-C07"], "Neurons (other)", STATE["C9"], r"SYNAP", "Down",
            lambda c: c == "C9", "synapse terms in c9")
add_pathway("C1_TYROBP_network_down", ["wang2024-C14"], "Microglia", "c1", r"TYROBP", "Down",
            lambda c: c == "C1", "WikiPathways TYROBP causal network in microglia")

# ---------------------------------------------------------------- named genes (text)
add_text("DA_TF_named", ["wang2024-C06"], "named_genes", "DA neurons", "", "", ["MYT1L", "SOX6"])
add_text("RIT2_IHC_named", ["wang2024-C02"], "named_genes", "Neurons (other)", "RIT2+ neurons", "", ["RIT2", "TH"])
add_text("metallothionein_named", ["wang2024-C09"], "named_genes", "All", "", "PD vs control",
         ["MT2A", "MT1E", "MT3"], "text list (up in PD)")
add_text("HSP_named", ["wang2024-C10"], "named_genes", "All", "", "PD vs control",
         ["HSPB1", "HSPH1", "HSPA1A", "HSPA1B", "HSP90AA1", "CRYAB"], "text list (up in PD)")
add_text("C7_3_DA_named", ["wang2024-C12"], "named_genes", "DA neurons", "c7_3", "PD vs control",
         ["SLC18A2", "ALDH1A1", "SLC6A3", "TH"], "text list (down in PD)")
add_text("shared_c9_c7_3_c6_2_named", ["wang2024-C13"], "named_genes", "DA neurons", "c9, c7_3, c6_2",
         "PD vs control", ["SSBP3", "PLCL1", "MTRNR2L1"], "text list (SSBP3 down; PLCL1, MTRNR2L1 up)")
add_text("PD_gene_expression_named", ["wang2024-C19"], "named_genes", "Microglia", "", "",
         ["LRRK2", "PRKN", "VPS13C", "DNAJC13"])
add_text("DA_familialPD_down_named", ["wang2024-C20"], "named_genes", "DA neurons", "c7_3 / c9", "PD vs control",
         ["SNCA", "DNAJC6", "CHCHD2"], "text list (down in PD)")
add_text("C9_familialPD_up_named", ["wang2024-C21"], "named_genes", "Neurons (other)", STATE["C9"], "PD vs control",
         ["DNAJC6", "UCHL1", "PRKN"], "text list (up in PD)")
add_text("SNCA_glia_named", ["wang2024-C22"], "named_genes", "Microglia", "", "PD vs control", ["SNCA"],
         "text list (up in PD microglia and oligodendrocytes)")
add_text("GWAS_neuron_named", ["wang2024-C24"], "named_genes", "Neurons (other)", "", "PD vs control",
         ["SV2C", "KTN1"], "text list (SV2C down in c7_3)")
add_text("C9_RIT2_named", ["wang2024-C25"], "named_genes", "Neurons (other)", STATE["C9"], "PD vs control",
         ["RIT2"], "text list (down in PD)")

# ---------------------------------------------------------------- unresolved
unres("wang2024-C15", "bulk-tissue PD meta-analysis DEGs (ref. 30)",
      "external bulk DEG list not in supplement; per-cluster DE signatures given instead")
unres("wang2024-C14", "IPA TYROBP causal network",
      "IPA network is proprietary; MSigDB WP_TYROBP_CAUSAL_NETWORK overlap genes from Table S6 used")
unres("wang2024-C24", "52 GWAS-locus DEGs (Fig. 6D)", "list only in figure, not in supplement tables")


# ---------------------------------------------------------------- write + checks
def write_and_check():
    OUT.mkdir(parents=True, exist_ok=True)
    sig = pd.DataFrame(rows, columns=COLS)
    un = pd.DataFrame(unresolved, columns=["claim_id", "term", "reason"])
    assert sig.gene.str.len().gt(0).all() and not sig.gene.isin(["nan", "None"]).any(), "empty gene"
    assert (sig.groupby("signature_id")["rank"].min() == 1).all(), "rank must start at 1"
    assert not sig.duplicated(["signature_id", "gene"]).any(), "duplicate gene in a signature"
    sig.to_csv(OUT / f"{PAPER}.tsv", sep="\t", index=False)
    un.to_csv(OUT / f"{PAPER}_unresolved.tsv", sep="\t", index=False)

    print(sig.groupby(["signature_id", "kind"]).size().rename("n_genes").to_string())
    print("\nsignatures by kind:", sig.drop_duplicates("signature_id").kind.value_counts().to_dict())
    datelike = sig.gene[sig.gene.str.match(r"^\d{1,2}-[A-Za-z]{3}$|^[A-Za-z]{3}-\d{1,2}$|^\d{4}-\d{2}-\d{2}")]
    print("date-like gene symbols:", sorted(set(datelike)) or "none")

    claims = yaml.safe_load(open(ROOT / "claims" / f"{PAPER}.yaml"))
    need = [c["id"] for c in claims if c["test_template"] in ("state_signature", "geneset_shift", "de_pseudobulk")
            or (c["test_template"].startswith("composition") and (c.get("cell_state") or "").strip())]
    in_sig = set(";".join(sig.claim_ids).split(";"))
    in_un = set(un.claim_id)
    missing = [c for c in need if c not in in_sig and c not in in_un]
    print(f"\ncoverage: {len(need)} required claims; {len([c for c in need if c in in_sig])} in signatures; "
          f"{len([c for c in need if c not in in_sig and c in in_un])} unresolved only; missing: {missing or 'none'}")
    assert not missing


write_and_check()

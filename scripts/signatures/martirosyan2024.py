"""Gene signatures for martirosyan2024 claims (see claims/SIGNATURES.md).

Sources (raw/supplements/martirosyan2024/):
  MOESM5  Supplementary Table 5: markers of major cell types and of subpopulations
          (Wilcoxon, Bonferroni-adjusted p, ln fold change)
  MOESM7  Supplementary Table 7: PD-vs-control DEGs per subpopulation (LMM with donor RE,
          BH-adjusted p < 0.1 as reported by the authors)
  MOESM10 Supplementary Table 10: 28 marker genes shared by the TH-enriched subpopulations

Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/martirosyan2024.py
"""
from pathlib import Path

import pandas as pd
import yaml

PAPER = "martirosyan2024"
ROOT = Path(__file__).resolve().parents[2]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

SUPP = cfg.RAW / "supplements" / PAPER
F_MK = SUPP / "13024_2023_699_MOESM5_ESM.xlsx"
F_DE = SUPP / "13024_2023_699_MOESM7_ESM.xlsx"
F_TH = SUPP / "13024_2023_699_MOESM10_ESM.xlsx"
OUT = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison",
        "gene", "rank", "stat", "padj", "source", "selection"]

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


C = lambda *n: [f"{PAPER}-C{i:02d}" for i in n]
BROAD = {"neurons": "Neurons (other)", "astrocytes": "Astrocytes", "microglia": "Microglia",
         "oligos": "Oligodendrocytes"}
# subpopulations whose broad class is more specific than the family
SPECIFIC = {"neurons0": "DA neurons", "neurons3": "Inhibitory neurons"}
STATE_NAME = {"neurons0": "Neurons0", "neurons2": "Neurons2", "neurons3": "Neurons3", "neurons4": "Neurons4",
              "astrocytes1": "Astrocytes1", "astrocytes2": "Astrocytes2 (TH-enriched)", "astrocytes3": "Astrocytes3",
              "astrocytes4": "Astrocytes4", "astrocytes5": "Astrocytes5", "microglia1": "Microglia1 (TH-enriched)",
              "microglia3": "Microglia3", "oligos0": "Oligos0", "oligos1": "Oligos1",
              "oligos2": "Oligos2 (TH-enriched)", "oligos3": "Oligos3", "oligos5": "Oligos5"}


def cell_type_of(sub):
    return SPECIFIC.get(sub, BROAD[sub.rstrip("0123456789")])


def read(path, sheet):
    d = pd.read_excel(path, sheet_name=sheet, header=2).dropna(how="all")
    return d


# ---------------------------------------------------------------- major cell-type markers (S5)
MK_COLS = ["cluster", "gene", "p", "lfc", "pct1", "pct2", "padj"]
SEL_MK = "padj<0.05 (Bonferroni) & ln-FC>0 vs other clusters, top 100 by ln-FC"
major = read(F_MK, "Major cell type markers (SN)")
major.columns = MK_COLS
major_claims = {"Astro": ("Astrocytes", C(2)), "Micro": ("Microglia", C(2)), "Oligo": ("Oligodendrocytes", C(2)),
                "OPC": ("OPCs", C(2)), "Neurons": ("Neurons (other)", C(1)), "T cells": ("T cells", C(3))}
for cl, (ct, cids) in major_claims.items():
    d = major[(major.cluster == cl) & (major.padj < 0.05) & (major.lfc > 0)].sort_values("lfc", ascending=False).head(100)
    add_table(f"major_{cl.replace(' ', '_')}_markers", cids, "state_markers", ct, "glia" if cids == C(2) else "",
              "", d, "gene", "lfc", "padj", f"{F_MK.name}|Major cell type markers (SN)", SEL_MK)

# ---------------------------------------------------------------- subpopulation markers (S5)
sub_claims = {"neurons0": C(4, 8, 9, 28), "neurons3": C(5), "neurons2": C(6), "neurons4": C(7),
              "astrocytes2": C(11), "astrocytes3": C(12), "astrocytes4": C(13), "astrocytes1": C(14),
              "astrocytes5": C(15), "microglia1": C(16), "oligos2": C(19), "oligos5": C(20), "oligos1": C(21),
              "oligos3": C(22), "oligos0": C(23)}
for sheet in ["Neuron markers (SN)", "Astrocyte markers (SN)", "Microglia markers (SN)",
              "Oligodendrocyte markers (SN)"]:
    m = read(F_MK, sheet)
    m.columns = MK_COLS
    for sub, d in m.groupby("cluster", sort=True):
        if sub not in sub_claims:
            continue
        d = d[(d.padj < 0.05) & (d.lfc > 0)].sort_values("lfc", ascending=False).head(100)
        add_table(f"{sub}_markers", sub_claims[sub], "state_markers", cell_type_of(sub), STATE_NAME[sub], "",
                  d, "gene", "lfc", "padj", f"{F_MK.name}|{sheet}", SEL_MK)

# ---------------------------------------------------------------- shared TH-enriched markers (S10)
th = pd.read_excel(F_TH, sheet_name="Shared genes between TH+ cells", header=None)
th_genes = [g for g in th.iloc[2:, 0].dropna().astype(str).str.strip() if g]
add_table("TH_enriched_shared_markers", C(25, 26), "state_markers", "All", "TH-enriched subpopulations", "",
          pd.DataFrame({"gene": th_genes}), "gene", None, None,
          f"{F_TH.name}|Shared genes between TH+ cells", "authors' list of shared markers (order as given)")

# ---------------------------------------------------------------- PD-vs-control DEGs (S7)
de_claims = {"neurons0": C(8, 9), "neurons3": C(10), "microglia1": C(17), "microglia3": C(18), "oligos2": C(24)}
for sheet in ["DEGs in neurons", "DEGs in astrocytes", "DEGs in microglia", "DEGs in oligodendrocytes"]:
    de = read(F_DE, sheet)
    de.columns = ["gene", "coef", "p", "z", "padj", "sub"]
    for sub, d in de.groupby("sub", sort=True):
        for kind, sign in (("de_up", 1), ("de_down", -1)):
            dd = d[(d.padj < 0.1) & (sign * d.coef > 0)].sort_values("coef", ascending=(sign < 0))
            if dd.empty:
                continue
            add_table(f"{sub}_PDvsCtrl_{kind[3:]}", de_claims.get(sub, []), kind, cell_type_of(sub),
                      STATE_NAME.get(sub, sub.capitalize()), "PD vs control", dd, "gene", "coef", "padj",
                      f"{F_DE.name}|{sheet}", "all DEGs in Table S7 (LMM with donor RE, BH padj<0.1), ranked by coef")

# ---------------------------------------------------------------- named genes (text)
add_text("neurons0_sugar_proteasome_named", C(8), "named_genes", "DA neurons", "Neurons0", "PD vs control",
         ["OGA", "CHST1", "KBTBD6"], "text list (down in PD)")
add_text("neurons0_golgi_named", C(9), "named_genes", "DA neurons", "Neurons0", "PD vs control",
         ["RAB6B", "RAB8B", "COG4"], "text list (dysregulated in PD)")
add_text("neurons3_named", C(10), "named_genes", "Inhibitory neurons", "Neurons3", "PD vs control",
         ["NEAT1", "CA2"], "text list (up in PD)")
add_text("microglia1_acetyltransferase_named", C(17), "named_genes", "Microglia", "Microglia1 (TH-enriched)",
         "PD vs control", ["KAT2B", "ATF2", "MCM3AP"], "text list (up in PD)")
add_text("microglia3_GTPase_named", C(18), "named_genes", "Microglia", "Microglia3", "PD vs control",
         ["DOCK1", "AGFG1", "ASAP1"], "text list (down in PD)")
add_text("oligos2_named", C(24), "named_genes", "Oligodendrocytes", "Oligos2 (TH-enriched)", "PD vs control",
         ["RBM22", "PRPF8", "HNRNPU", "WNK1", "ENSA", "CDKN1B", "CDKN1C", "BIRC2", "SYVN1", "CAMK2G"],
         "text list (up in PD)")
add_text("TH_in_situ_named", C(27), "named_genes", "All", "TH+ glia", "", ["TH", "SLC6A3", "SNAP25"])
add_text("TH_glia_Kamath_named", C(26), "named_genes", "All", "TH+ glia", "PD vs control", ["TH"])
add_text("monogenic_PD_neurons0_named", C(28), "named_genes", "DA neurons", "Neurons0", "",
         ["PARK7", "PINK1", "ATP13A2", "VPS35", "SYNJ1", "SNCA", "DNAJC6"])
add_text("monogenic_PD_microglia_named", C(29), "named_genes", "Microglia", "", "",
         ["LRRK2", "VPS13C", "SNCA", "DNAJC6"])

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
    datelike = sig[sig.gene.str.match(r"^\d{1,2}-[A-Za-z]{3}$|^[A-Za-z]{3}-\d{1,2}$|^\d{4}-\d{2}-\d{2}|^\d{5}$")]
    print("date-like gene symbols (kept as reported):",
          datelike[["signature_id", "gene"]].values.tolist() or "none")

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

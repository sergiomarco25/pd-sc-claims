"""Genetic enrichment claims of the PD field: PD GWAS (Nalls 2019) MAGMA gene-set tests on each dataset's gene sets
(scripts/16a-16c), and familial PD genes in each class's specific genes. The engine's claim re-tests call
`genetic_test` for every genetic_enrichment claim ([tests] genetic_test in config/project.toml).
"""
import sys
from functools import lru_cache
from pathlib import Path

import pandas as pd
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from claimstore.retest import RES, summarise  # noqa: E402

GWAS_SETS = RES / "gwas" / "sets"


@lru_cache(maxsize=32)
def gsa(dataset):
    f = GWAS_SETS / f"{dataset}.gsa.out"
    if not f.exists():
        return None
    t = pd.read_csv(f, sep=r"\s+", comment="#")
    name = t["FULL_NAME"] if "FULL_NAME" in t else t["VARIABLE"]
    return t.assign(set=name.str.replace("_", " ")).set_index("set")


def genetic_test(ctx, dataset, statement, modality, cell_state):
    """PD GWAS (Nalls 2019) MAGMA competitive gene-set tests on this dataset's sets (scripts/16c_gwas_sets.py)."""
    st = statement.lower()
    na = lambda why: [{"template_used": "gwas_magma", "outcome": "not_testable", "note": why}]
    if any(w in st for w in ["alzheimer", "diabetes", "t2d", "age of onset", "east asian"]):
        return na("needs a GWAS other than PD (Nalls 2019)")
    if any(w in st for w in ["atac", "cre", "peak", "accessib", "chromatin"]) or "atac" in str(modality).lower():
        return na("claim is about chromatin regions (needs ATAC)")
    if str(cell_state).strip() and "familial" not in st:
        return na("subtype gene set is fixed by the source paper (verified once, not replicable per dataset)")
    classes = [c for c in str(ctx.cell_classes).split(";") if c and c != "ALL"]
    if "familial" in st:
        f = GWAS_SETS / f"{dataset}.familial.tsv"
        if not f.exists():
            return na("no familial-gene overlap file")
        t = pd.read_csv(f, sep="\t").set_index("cell_class")
        rows = [{"region": "all", "cell_class": c, "z": stats.norm.isf(max(t.loc[c, "p"], 1e-300)), "p": t.loc[c, "p"],
                 "effect": t.loc[c, "familial_in_set"], "n_case": int(t.loc[c, "n_genes"]), "n_ref": 0}
                for c in classes if c in t.index]
        return [summarise(pd.DataFrame(rows), 1, "familial_overlap")]
    t = gsa(dataset)
    if t is None:
        return na("no MAGMA gene-set results for this dataset")
    if any(w in st for w in ["upregulated", "deg", "differentially"]):
        con = "PDD_vs_control" if "pdd" in st else "LRRK2_vs_control" if "lrrk2" in st else "PDany_vs_control"
        keys = [f"deup::{con}::{c}" for c in classes]
        label = "gwas_magma_de_up"
    else:
        keys = [f"spec::{c}" for c in classes]
        label = "gwas_magma_specificity"
    null_claim = any(w in st for w in ["no significant", "not significantly", "no enrichment"])
    if null_claim:  # "no cell type is enriched": reproduced if no specificity set passes Bonferroni
        sp = t[t.index.str.startswith("spec::")]
        if sp.empty:
            return na("no specificity sets")
        pmin = sp.P.min() * len(sp)
        return [{"template_used": "gwas_magma_specificity", "outcome": "reproduced" if pmin >= 0.05 else "contradicted",
                 "z": float(stats.norm.isf(min(sp.P.min(), 1))), "p": float(min(pmin, 1)), "n_tests": len(sp),
                 "classes_tested": ";".join(sp.index), "note": "null claim: no class enriched after Bonferroni"}]
    rows = [{"region": "all", "cell_class": k.split("::")[-1], "z": stats.norm.isf(max(t.loc[k, "P"], 1e-300)),
             "p": float(t.loc[k, "P"]), "effect": float(t.loc[k, "BETA"]), "n_case": int(t.loc[k, "NGENES"]), "n_ref": 0}
            for k in keys if k in t.index]
    return [summarise(pd.DataFrame(rows), 1, label)]

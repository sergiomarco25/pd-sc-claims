"""Evidence map: the 456 claims grouped into biological themes, with how each theme held up across datasets.

Usage: python -I scripts/26_evidence_map.py
Theme assignment is rule-based (regular expressions on the claim statement plus its cell type; a claim can belong to
several themes). Writes results/biology/claim_themes.tsv (claim x theme) and evidence_map.tsv (per theme: claims,
papers, robustness counts, share robust or supported, internal reproduction).
"""
from pathlib import Path
import re
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[1]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

OUT = cfg.RESULTS / "biology"
THEMES = {  # theme: (statement pattern, cell types that also qualify)
    "DA neuron vulnerability": (r"dopamin|\bDA\b|DaN|AGTR1|SOX6|CALB1|\bTH\b|neuromelanin|RIT2", {"DA neurons"}),
    "Microglial activation": (r"microgli|GPNMB|DAM\b|IL1B|CD74|macrophage|myeloid", {"Microglia"}),
    "Astrocyte reactivity": (r"astrocyt|GFAP|CD44|reactive astro|AR1|AR2", {"Astrocytes"}),
    "Oligodendroglia & myelin": (r"oligodend|\bOPC|myelin|PDAO|OPALIN|ODC", {"Oligodendrocytes", "OPCs"}),
    "Proteostasis / heat shock / UPR": (r"heat.?shock|\bHSP|HSPA|DNAJ|unfolded|\bUPR\b|chaperon|proteostas|protein folding|CRYAB", set()),
    "Mitochondria & metabolism": (r"mitochond|OXPHOS|oxidative phosph|metaboli|glycoly|lipid|cholesterol|insulin", set()),
    "Adaptive immunity / T cells": (r"T cell|T-cell|lymphocyt|CD8|CD4|TCR|clonal|NK cell", {"T cells"}),
    "Vasculature & barrier": (r"endothel|pericyte|vascular|fibroblast|VLMC|blood.?brain|coagulation", {"Endothelial", "Pericytes", "Fibroblasts"}),
    "Synaptic & neuronal function": (r"synap|excitatory|inhibitory|interneuron|MSN|glutamat|GABA", {"Excitatory neurons", "Inhibitory neurons", "MSNs"}),
    "Genetic risk convergence": (r"GWAS|heritab|familial|risk (gene|variant|loci)|monogenic|LRRK2|GBA", set()),
    "Lewy pathology & stage": (r"Braak|Lewy|stage|pathology|alpha-?synuclein|SNCA", set()),
}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    claims = pd.DataFrame(yaml.safe_load((ROOT / "claims" / "claims.yaml").read_text())).set_index("id")
    r4 = pd.read_csv(ROOT / "reports" / "R4_cross_claims.tsv", sep="\t").set_index("claim_id")
    r3 = pd.read_csv(ROOT / "reports" / "R3_internal_claims.tsv", sep="\t").set_index("claim_id")
    rows = []
    for cid, c in claims.iterrows():
        cts = set(c.cell_type if isinstance(c.cell_type, list) else [c.cell_type])
        for th, (pat, cells) in THEMES.items():
            if re.search(pat, str(c.statement), re.I) or cts & cells:
                rows.append({"claim_id": cid, "theme": th, "paper": cid.rsplit("-C", 1)[0], "direction": c.direction,
                             "robustness": r4.robustness.get(cid, "insufficient"), "internal": r3.outcome.get(cid, "not_testable"),
                             "statement": c.statement})
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "claim_themes.tsv", sep="\t", index=False)
    tested = t[t.robustness != "insufficient"]
    g = t.groupby("theme")
    em = pd.DataFrame({"claims": g.size(), "papers": g.paper.nunique()})
    em = em.join(pd.crosstab(tested.theme, tested.robustness)).fillna(0)
    em["tested"] = tested.groupby("theme").size()
    em["share_robust_or_supported"] = (em.get("robust", 0) + em.get("supported", 0)) / em.tested
    em["share_internal_reproduced"] = t[t.internal != "not_testable"].groupby("theme").internal.apply(lambda s: (s == "reproduced").mean())
    em = em.sort_values("share_robust_or_supported", ascending=False)
    em.to_csv(OUT / "evidence_map.tsv", sep="\t")
    print(em.round(2).to_string())


if __name__ == "__main__":
    main()

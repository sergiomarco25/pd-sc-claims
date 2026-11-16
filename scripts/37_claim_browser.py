"""Interactive claim browser: every claim ranked by cross-dataset evidence, with region and donor context.

Usage: python -I scripts/37_claim_browser.py
Inputs: claims/claim_evidence.tsv (scripts/33) and results/meta/claim_dataset_matrix.tsv.
Writes reports/claim_browser.html (self-contained; data embedded as JSON; template reports/claim_browser_template.html).
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

MAT = cfg.RESULTS / "meta/claim_dataset_matrix.tsv"
DS = ["smajic2022", "kamath2022", "lee2023", "martirosyan2024", "wang2024", "adams2024", "ma2025", "volpato2025",
      "menon2026", "barbareyes2026", "xu2023", "prashant2024", "feleke2021", "zhu2024", "chia2025", "jin2024",
      "glock2026", "lee2026"]
PAPER_DS = {"barbareyes2025": "barbareyes2026", "clarence2025": "prashant2024"}
CODE = {"reproduced": "R", "direction_only": "D", "not_reproduced": "N", "contradicted": "C"}
REGION_GROUP = {"SN": "midbrain / SN", "midbrain": "midbrain / SN", "putamen": "striatum", "caudate": "striatum",
                "striatum": "striatum", "DLPFC": "cortex", "PFC": "cortex", "cingulate": "cortex", "temporal": "cortex",
                "multiple": "several regions", "cerebellum": "other"}


def clean(v):
    if isinstance(v, float) and not np.isfinite(v):
        return None
    if isinstance(v, (np.integer,)):
        return int(v)
    if isinstance(v, (np.floating,)):
        return float(v)
    return v


def main():
    e = pd.read_csv(ROOT / "claims" / "claim_evidence.tsv", sep="\t")
    m = pd.read_csv(MAT, sep="\t").drop_duplicates("claim_id").set_index("claim_id")
    recs = []
    for r in e.itertuples():
        src = PAPER_DS.get(r.paper, r.paper)
        mm = m.loc[r.claim_id] if r.claim_id in m.index else pd.Series(dtype=object)
        # one letter per dataset; lower case marks the source dataset (own data), "_" an untested source
        out = "".join(CODE.get(r.own_data_outcome, "_").lower() if ds == src else CODE.get(mm.get(ds, ""), "-")
                      for ds in DS)
        rec = {k: clean(getattr(r, k)) for k in
               ["rank", "tier", "kind", "claim_id", "paper", "year", "statement", "themes", "claim_type", "direction",
                "comparison", "region", "cell_type", "unit_of_analysis", "reported_stat", "source_regions",
                "source_cases", "source_controls", "own_data_outcome", "k", "n_reproduced", "n_direction_only",
                "n_not_reproduced", "n_contradicted", "z", "p", "I2", "same_region_k", "same_region_rep",
                "other_region_k", "other_region_rep", "donors_cases_tested", "donors_controls_tested",
                "regions_tested", "datasets_reproducing", "datasets_contradicting", "quote", "location"]}
        rec["outcomes"] = out
        rec["rgroup"] = REGION_GROUP.get(r.region, "other")
        recs.append(rec)
    tpl = (ROOT / "reports" / "claim_browser_template.html").read_text()
    html = tpl.replace("/*__DATA__*/[]", json.dumps(recs, ensure_ascii=False)).replace(
        "/*__DS__*/[]", json.dumps(DS))
    (ROOT / "reports" / "claim_browser.html").write_text(html)
    print(len(recs), "claims;", round(len(html) / 1e6, 2), "MB")


if __name__ == "__main__":
    main()

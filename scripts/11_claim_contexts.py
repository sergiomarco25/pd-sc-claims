"""Map every claim to the data needed to test it: dataset, regions, cell classes and contrast.

Usage: python -I scripts/11_claim_contexts.py
Writes claims/contexts.tsv (one row per claim; the source dataset for internal replication) and prints claims whose
comparison text matched no rule. Rules are applied in order; claims/contrast_overrides.tsv (claim_id, contrast,
note) takes precedence over the rules.

Contrasts (case vs reference, in pdrep.harmonise disease labels):
  PDany_vs_control   PD+PDD vs control (papers that say "PD" usually pooled PD with dementia)
  LBDany_vs_control  PD+PDD+DLB+LBD vs control
  PD_vs_control, PDD_vs_control, DLB_vs_control, AD_vs_control, iLBD_vs_control
  PDD_vs_PD, DLB_vs_PD, DLB_vs_PDD
  PD_vs_agedcontrol, aged_vs_youngcontrol            (adams2024 age_group)
  GBA_vs_control, LRRK2_vs_control                  (dehestani2024 genotype)
  LBearly/LBmid/LBlate_vs_control   Lewy-body stage groups (Braak 1-2 / 3-4 / >=5) vs controls
  braak            association with Lewy-body stage across donors
  sex              male vs female
  within           no disease contrast (state existence, markers, cell-type specificity)
  gliad_groups     GLIAD donor groups (UPR-high / Reactive vs homeostatic), defined per dataset by the GLIAD scores
  PD_AD_concordance  agreement of PD-vs-control and AD-vs-control effects (datasets with AD donors)
  not_available    the comparison needs data no dataset has (spatial tier, white matter, haplotype)
"""
from pathlib import Path
import re
import pandas as pd
import yaml
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pdrep.regions import REGION_MAP  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]

from pdrep import config as cfg  # noqa: E402

H = cfg.HARMONISED
PAPER_TO_DATASET = {"barbareyes2025": "barbareyes2026", "clarence2025": "prashant2024", "fullard2025": "lee2026"}
RULES = [  # (regex on lower-cased comparison, contrast)
    (r"aged control vs young|aged vs young|young control|cpp", "aged_vs_youngcontrol"),
    (r"pd vs aged control", "PD_vs_agedcontrol"),
    (r"gba1?-pd vs control|gba1-pd", "GBA_vs_control"),
    (r"lrrk2-pd", "LRRK2_vs_control"),
    (r"dlb vs pdd", "DLB_vs_PDD"),
    (r"pdd vs pd\b", "PDD_vs_PD"),
    (r"dlb vs pd\b", "DLB_vs_PD"),
    (r"pd \(braak 5-6\) vs control|braak 4-6 vs control|braak lb 3-6 vs control", "LBlate_vs_control"),
    (r"braak lb 1-2 vs control", "LBearly_vs_control"),
    (r"braak|lb stage|lewy body pathology score|lewy body \+ neurite|pseudotime", "braak"),
    (r"male vs female", "sex"),
    (r"^pdd vs control|^pdd/dlbd vs control|^pdd, dlb", "PDD_vs_control"),
    (r"^dlb vs control", "DLB_vs_control"),
    (r"lbd vs control|pd/lbd|pd/ilbd", "LBDany_vs_control"),
    (r"^pd vs control|^pd effect|^pd, pdd, dlb vs control", "PDany_vs_control"),
    (r"^ad vs control", "AD_vs_control"),
    (r"within|across|vs other|vs reference|vs mol-|vs opcs-|vs dans|vs healthy odc|specificity|markers|"
     r"all groups|conserved|subtypes|clusters|cell types|neighbourhood|^na$|presence", "within"),
]
# "No difference" claims (scored with an equivalence-style rule); GWAS null claims are handled in the GWAS test.
NULL_PAT = (r"not significantly|no significant|not changed|unchanged|no difference|not (?:reduced|increased|altered|"
            r"different|depleted|enriched)|comparable|did not differ|does not differ|no change")
NEURONS = ["DA neurons", "MSNs", "Excitatory neurons", "Inhibitory neurons", "Neurons (other)"]


def contrast_of(text):
    t = str(text).strip().lower()
    for pat, con in RULES:
        if re.search(pat, t):
            return con
    return ""


def as_list(v):
    return v if isinstance(v, list) else ([] if v in (None, "") else [v])


def main():
    claims = pd.DataFrame(yaml.safe_load((ROOT / "claims" / "claims.yaml").read_text())).fillna("")
    ov_path = ROOT / "claims" / "contrast_overrides.tsv"
    ov = pd.read_csv(ov_path, sep="\t").set_index("claim_id").contrast.to_dict() if ov_path.exists() else {}
    regions_by_ds = {f.name[:-len("_donors.tsv")]: set(pd.read_csv(f, sep="\t").region.astype(str))
                     for f in H.glob("*_donors.tsv")}
    rows = []
    for c in claims.itertuples():
        ds = PAPER_TO_DATASET.get(c.paper, c.paper)
        con = ov.get(c.id, contrast_of(c.comparison))
        have = regions_by_ds.get(ds, set())
        want = set().union(*[REGION_MAP.get(r, have) for r in as_list(c.region)]) if as_list(c.region) else have
        regs = sorted(have & want) or sorted(have)  # claim region absent from its own data -> use what exists
        cts = as_list(c.cell_type)
        classes = sorted({x for ct in cts for x in (NEURONS if ct == "Neurons (other)" else ["ALL"] if ct == "All" else [ct])})
        status = ("no dataset" if ds not in regions_by_ds else "not_testable" if c.test_template == "not_testable"
                  else "not_available" if con == "not_available" else "no contrast" if not con else "ok")
        rows.append({"claim_id": c.id, "paper": c.paper, "dataset": ds, "regions": ";".join(regs),
                     "cell_classes": ";".join(classes), "cell_state": c.cell_state, "contrast": con,
                     "direction": c.direction, "claim_type": c.claim_type, "test_template": c.test_template,
                     "status": status, "comparison_text": c.comparison,
                     "null_claim": bool(str(c.direction) in ("NA", "", "None") and c.test_template != "genetic_enrichment"
                                        and re.search(NULL_PAT, str(c.statement), re.I))})
    ctx = pd.DataFrame(rows)
    ctx.to_csv(ROOT / "claims" / "contexts.tsv", sep="\t", index=False)
    print(ctx.status.value_counts().to_string(), "\n")
    print(ctx[ctx.status == "ok"].contrast.value_counts().to_string(), "\n")
    un = ctx[ctx.status == "no contrast"]
    print(f"{len(un)} claims with no contrast rule:")
    print(un[["claim_id", "comparison_text"]].to_string(index=False, max_colwidth=90))


if __name__ == "__main__":
    main()

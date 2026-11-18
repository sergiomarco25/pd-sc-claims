"""Map every claim to the data needed to test it: dataset, regions, cell classes and contrast.

Run through the project's entry point (e.g. python -I scripts/11_claim_contexts.py).
Writes claims/contexts.tsv (one row per claim; the source dataset for internal replication) and prints claims whose
comparison text matched no rule. The rules, the no-difference pattern, how a claim cell type expands to tested cell
classes and how claim regions map to dataset regions come from config/contexts.toml; claims/contrast_overrides.tsv
(claim_id, contrast, note) takes precedence over the rules. Contrast codes are those of config/contrasts.toml plus
the non-group codes listed in [contrast_labels] (e.g. a stage covariate, `within`, `not_available`).
"""
import re
import tomllib

import pandas as pd
import yaml

from claimstore import config as cfg

ROOT = cfg.REPO
H = cfg.HARMONISED
PAPER_TO_DATASET = cfg.PAPER_TO_DATASET  # [claims] paper_to_dataset in config/project.toml
SPEC_FILE = cfg.REPO / "config" / "contexts.toml"
SPEC = tomllib.loads(SPEC_FILE.read_text()) if SPEC_FILE.is_file() else {}
RULES = [tuple(r) for r in SPEC.get("rules", [])]  # (regex on lower-cased comparison, contrast)
# "No difference" claims (scored with an equivalence-style rule)
NULL_PAT = SPEC.get("null_pattern", r"(?!x)x")
EXPAND = {k: list(v) for k, v in SPEC.get("cell_class_expand", {}).items()}
REGION_MAP = {k: set(v) for k, v in SPEC.get("region_map", {}).items()}
CONTRAST_LABELS = {k: tuple(v) for k, v in SPEC.get("contrast_labels", {}).items()}


def target_regions(claim_regions, have):
    """Regions of a dataset to test a claim in, and whether they match the claim's region ('same' / 'other')."""
    claim_regions = [r for r in claim_regions if r]
    if not claim_regions or "multiple" in claim_regions:
        return sorted(have), "same"
    want = set().union(*[REGION_MAP.get(r, set()) for r in claim_regions])
    hit = sorted(have & want)
    return (hit, "same") if hit else (sorted(have), "other")


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
        classes = sorted({x for ct in cts for x in EXPAND.get(ct, [ct])})
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

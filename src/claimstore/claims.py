"""Validate per-paper claim files and merge them into claims/claims.yaml and claims/claims.tsv.

Run through the project's entry point (e.g. python -I scripts/05_merge_claims.py).
Checks every claim against the project's claim schema (config/schema.toml): required fields, controlled
vocabularies, id format, unique ids.
Problems are listed in claims/validation.tsv (claim id, field, problem); values are not changed here.
"""
from pathlib import Path
import re
import pandas as pd
import yaml

from claimstore import config as cfg
from claimstore import schema  # config/schema.toml: the same vocabularies as the JSON Schema

CL = cfg.REPO / "claims"
REQUIRED = schema.REQUIRED
# fields checked here against their vocabulary (regions and cell types are checked by the JSON Schema and the contexts)
VOCAB = {k: set(schema.VOCAB[k]) for k in ("centrality", "claim_type", "direction", "testable_with", "test_template",
                                           "modality")}
LISTY = {"genes", "region", "cell_type", "modality"}


def as_list(v):
    if v is None:
        return []
    return v if isinstance(v, list) else [x.strip() for x in str(v).replace(";", ",").split(",") if x.strip()]


def main():
    claims, problems = [], []
    for f in sorted(CL.glob("*.yaml")):
        if f.stem in ("claims", "findings", "background", "questions"):  # not per-paper claim files
            continue
        pid = f.stem
        data = yaml.safe_load(f.read_text()) or []
        if not isinstance(data, list):
            problems.append((pid, "file", "not a YAML list")); continue
        for c in data:
            cid = c.get("id", f"{pid}-?")
            for k in REQUIRED:
                if c.get(k) in (None, "", []):
                    problems.append((cid, k, "missing"))
            if not re.fullmatch(rf"{re.escape(pid)}-C\d\d", str(cid)):
                problems.append((cid, "id", f"expected {pid}-Cnn"))
            if c.get("paper") != pid:
                problems.append((cid, "paper", f"{c.get('paper')} != file {pid}"))
            for k, allowed in VOCAB.items():
                for v in (as_list(c.get(k)) if k in LISTY else [c.get(k)]):
                    if v not in (None, "") and str(v) not in allowed:
                        problems.append((cid, k, f"'{v}' not in vocabulary"))
            claims.append(c)
    ids = [c.get("id") for c in claims]
    for d in {i for i in ids if ids.count(i) > 1}:
        problems.append((d, "id", "duplicate"))

    (CL / "claims.yaml").write_text(yaml.safe_dump(claims, sort_keys=False, allow_unicode=True, width=120))
    flat = pd.DataFrame([{k: "; ".join(map(str, v)) if isinstance(v, list) else v for k, v in c.items()} for c in claims])
    flat.to_csv(CL / "claims.tsv", sep="\t", index=False)
    pd.DataFrame(problems, columns=["claim", "field", "problem"]).to_csv(CL / "validation.tsv", sep="\t", index=False)
    print(f"{len(claims)} claims from {flat.paper.nunique() if len(flat) else 0} papers; {len(problems)} problems")
    if len(flat):
        print(flat.groupby("paper").size().to_string())
        for col in ["centrality", "claim_type", "test_template", "testable_with", "unit_of_analysis"]:
            print("\n" + flat[col].value_counts().to_string())


if __name__ == "__main__":
    main()

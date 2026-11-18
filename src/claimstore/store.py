"""Findings and questions as claims (claims/ROADMAP.md sections 2 and 4).

A *finding* is a claim made by our own analyses. It uses the claim schema of claims/SCHEMA.md plus:
  origin       "finding"
  analysis     {script, results, date}: what produced it
  evidence     a reference to a results table, resolved at build time (51b), so the evidence follows reruns:
                 table     path under results/ (e.g. biology/programme_effects.tsv)
                 filter    {column: value or [values]}: rows that belong to this finding
                 columns   mapping to dataset, est, se and n (or n_case and n_ref); defaults below
                 scale     smd (Hedges' g), fisher_z (correlation), logit (composition), or other
                 dataset   constant dataset id when the table has no dataset column (single-cohort results)
  contrast     a contrast code of the store (51b CONTRASTS), or a covariate code
  null_claim   true when the finding is "no difference / no association"

Analysis scripts call add_finding() at their end; findings live in claims/findings.yaml (git), keyed by id, so a
rerun updates the entry instead of duplicating it. Questions (claims/questions.yaml) follow the same pattern.
"""
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

ROOT = Path(__file__).resolve().parents[2]
FINDINGS = ROOT / "claims" / "findings.yaml"
QUESTIONS = ROOT / "claims" / "questions.yaml"

from claimstore import config as cfg  # noqa: E402

RESULTS = cfg.RESULTS
DEFAULT_COLS = {"dataset": "dataset", "est": "est", "se": "se", "n": "n", "n_case": "n_case", "n_ref": "n_ref",
                "region": "region"}
REQUIRED = ["id", "statement", "claim_type", "direction", "cell_type", "region", "contrast", "evidence", "analysis"]


def _load(path):
    return yaml.safe_load(path.read_text()) or [] if path.exists() else []


def _dump(path, items):
    path.write_text(yaml.safe_dump(items, sort_keys=False, allow_unicode=True, width=120))


def add_finding(id, statement, claim_type, direction, cell_type, region, contrast, evidence, analysis,
                cell_state="", modality="snRNA", comparison="", null_claim=False, centrality="main", genes=None,
                method="", notes="", themes=""):
    """Insert or update one finding in claims/findings.yaml (matched on id). Returns the stored entry."""
    entry = dict(id=id, origin="finding", paper="this-study", statement=statement, quote="",
                 location=f"{analysis.get('script', '')} -> {analysis.get('results', '')}", centrality=centrality,
                 modality=modality, region=region, cell_type=cell_type, cell_state=cell_state, claim_type=claim_type,
                 direction=direction, comparison=comparison or contrast, genes=genes or [], gene_source="",
                 reported_stat="", unit_of_analysis="donors", method=method or "donor-level meta-analysis (this study)",
                 testable_with="rna", test_template="finding", requires="", notes=notes, themes=themes,
                 contrast=contrast, null_claim=bool(null_claim), evidence=evidence,
                 analysis={**analysis, "date": analysis.get("date", str(date.today()))})
    missing = [k for k in REQUIRED if entry.get(k) in (None, "", [])]
    if missing:
        raise ValueError(f"finding {id}: missing {missing}")
    items = [x for x in _load(FINDINGS) if x.get("id") != id]
    items.append(entry)
    _dump(FINDINGS, sorted(items, key=lambda x: x["id"]))
    return entry


def add_question(id, text, asked_by="", status="open", answer="", records=None, proposed=None):
    """Insert or update a question (claims/questions.yaml). `proposed`: a finding draft (add_finding arguments)."""
    entry = dict(id=id, text=text, asked_by=asked_by, date=str(date.today()), status=status, answer=answer,
                 records=records or [], proposed=proposed or {})
    items = [x for x in _load(QUESTIONS) if x.get("id") != id]
    items.append(entry)
    _dump(QUESTIONS, sorted(items, key=lambda x: x["id"]))
    return entry


def _to_smd(est, se, scale, n_case, n_ref, n):
    """Effect and SE on the standardised-mean-difference scale used by the rating."""
    if scale == "smd":
        return est, se
    if scale == "fisher_z":
        r = np.tanh(est)
        d = 2 * r / np.sqrt(max(1 - r ** 2, 1e-6))
        return d, se * 2 / max(1 - r ** 2, 1e-6) ** 1.5 * (1 - r ** 2)  # delta method through r = tanh(z)
    z = est / se if se else np.nan  # other scales: via the test statistic and group sizes
    if n_case and n_ref:
        k = np.sqrt(1 / n_case + 1 / n_ref)
        d = z * k
        return d, np.sqrt(k ** 2 + d ** 2 / (2 * (n_case + n_ref)))
    if n:
        k = 2 / np.sqrt(n)
        return z * k, np.sqrt(k ** 2)
    return np.nan, np.nan


def resolve_evidence(ev):
    """Rows of a results table -> one test per dataset: dataset, est, se (native scale), d, se_d (SMD), n_case, n_ref,
    regions. Regions of one dataset are combined by fixed-effect inverse-variance weighting."""
    t = pd.read_csv(RESULTS / ev["table"], sep="\t", low_memory=False)
    for col, val in (ev.get("filter") or {}).items():
        t = t[t[col].isin(val if isinstance(val, list) else [val])]
    cols = {**DEFAULT_COLS, **(ev.get("columns") or {})}
    if ev.get("dataset"):
        t = t.assign(**{cols["dataset"]: ev["dataset"]})
    t = t[t[cols["dataset"]].notna() & t[cols["est"]].notna() & t[cols["se"]].notna()]
    rows = []
    for ds, g in t.groupby(cols["dataset"]):
        w = 1 / g[cols["se"]].astype(float) ** 2
        est = float((w * g[cols["est"]]).sum() / w.sum())
        se = float(np.sqrt(1 / w.sum()))
        nc = float(g[cols["n_case"]].max()) if cols["n_case"] in g else np.nan
        nr = float(g[cols["n_ref"]].max()) if cols["n_ref"] in g else np.nan
        n = float(g[cols["n"]].max()) if cols["n"] in g else (nc + nr if np.isfinite(nc + nr) else np.nan)
        if not np.isfinite(nc) and np.isfinite(n):  # tables with a total n only: carry it as the donor count
            nc, nr = n, np.nan
        d, se_d = _to_smd(est, se, ev.get("scale", "smd"), nc, nr, n)
        regions = ";".join(sorted(g[cols["region"]].astype(str).unique())) if cols["region"] in g else ""
        rows.append(dict(dataset=str(ds), est=est, se=se, d=d, se_d=se_d, n_case=nc, n_ref=nr, n=n, regions=regions,
                         k_units=len(g)))
    return pd.DataFrame(rows)


UP = {"up", "enriched", "associated_pos", "exists"}
DOWN = {"down", "depleted", "associated_neg"}


def finding_tests(f, ev):
    """One claim_dataset row per dataset of a finding: resolved evidence (resolve_evidence) -> outcome as claimed.
    `f` is the findings.yaml entry; the outcome uses its direction and null_claim, p from the SMD z (two-sided)."""
    from scipy.stats import norm
    exp = 1 if f["direction"] in UP else -1 if f["direction"] in DOWN else 0
    tests = []
    for _, e in ev.iterrows():
        z = e.d / e.se_d if e.se_d else np.nan
        p = float(2 * norm.sf(abs(z))) if np.isfinite(z) else np.nan
        if f.get("null_claim"):
            outcome = "reproduced" if p >= 0.05 else "contradicted"
        elif exp == 0:
            outcome = "reproduced" if p < 0.05 else "not_reproduced"
        else:
            s_ = np.sign(z) * exp
            outcome = ("reproduced" if s_ > 0 else "contradicted") if p < 0.05 else ("direction_only" if s_ > 0 else "not_reproduced")
        tests.append(dict(claim_id=f["id"], dataset=e.dataset, role="evidence", engine="finding", variant="both",
                          contrast=f["contrast"], region_match="same", template_used="finding", outcome=outcome,
                          effect=e.est, z=z, p=p, n_case=e.n_case, n_ref=e.n_ref, regions_tested=e.regions,
                          classes_tested=f["cell_type"], set_test=f["evidence"]["table"], note="",
                          d_native=e.d, se_native=e.se_d))
    return tests

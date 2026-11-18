"""Compile the claim store: claims with resolved context, datasets/cohorts, claim <-> dataset evidence, relations.

Inputs (git): claims/claims.yaml, claims/contexts.tsv, claims/vocab/*.tsv (51a_vocab.py), claims/structure_overrides.tsv,
  claims/consensus_statements.tsv, claims/source_discrepancies.tsv, claims/claim_evidence.tsv, claims/signatures_all.tsv.gz,
  inventory/papers.tsv, inventory/cohorts.yaml (curated cohort context: ancestry, brain bank, sampling).
Inputs (lustre): harmonised/<dataset>_donors.tsv, results/tests/{internal,internal_gwas,cross,cross_gwas}/*.tsv.
Outputs:
  claims/store/claim_structure.tsv        derived subject / predicate / contrast / context per claim (review this)
  claims/store/claims_contextualised.jsonl one self-contained record per claim (for AI retrieval)
  claims/store/validation.tsv             problems found while building
  results/store/claim_store.sqlite        all tables (see claims/STORE_DESIGN.md)

Inputs that a new project may not have yet (no vocabularies, no tests, no evidence summary, no signatures, no
consensus statements, no harmonised datasets) are treated as empty. Wording about the field comes from the project
config: [claims] sample_type, heritability_label, stage_available_column and staging_field; contrast labels from
config/contexts.toml [contrast_labels].

Run through the project's entry point (e.g. python scripts/51b_claim_store.py, in the analysis environment).
"""
import json
import re
import sqlite3
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
import yaml

from claimstore import config as cfg
from claimstore import contrasts as _contrasts
from claimstore.contexts import CONTRAST_LABELS

ROOT = cfg.REPO
H, RES = cfg.HARMONISED, cfg.RESULTS
VOCAB, OUT_GIT, OUT_DB = ROOT / "claims" / "vocab", ROOT / "claims" / "store", RES / "store"
SAMPLE = cfg.CLAIMS.get("sample_type", "tissue")                       # e.g. post-mortem brain
HERITABILITY = cfg.CLAIMS.get("heritability_label", "genetic risk")     # subject of genetic-enrichment claims
STAGE_COL = cfg.CLAIMS.get("stage_available_column", "stage_available")  # dataset column: any donor staged
STAGING = cfg.CLAIMS.get("staging_field", "staging")  # cohort-file field (and dataset column) describing the staging

# contrast code (claims/contexts.tsv) -> (kind, case group, reference group): config/contexts.toml [contrast_labels]
CONTRASTS = CONTRAST_LABELS
# (claim_type, direction) -> predicate; then direction-only fallback
PRED_TYPE = {
    ("composition", "up"): "increased_proportion_in", ("composition", "enriched"): "increased_proportion_in",
    ("composition", "down"): "decreased_proportion_in", ("composition", "depleted"): "decreased_proportion_in",
    ("de_gene", "up"): "increased_expression_in", ("de_gene", "down"): "decreased_expression_in",
    ("de_gene", "NA"): "altered_expression_in", ("de_pathway", "up"): "increased_expression_in",
    ("de_pathway", "down"): "decreased_expression_in", ("de_pathway", "NA"): "altered_expression_in",
    ("genetic_enrichment", "enriched"): "heritability_enriched_in",
    ("genetic_enrichment", "associated_pos"): "heritability_enriched_in",
    ("genetic_enrichment", "down"): "heritability_not_enriched_in",
    ("genetic_enrichment", "NA"): "heritability_tested_in",
    ("chromatin", "up"): "increased_accessibility_in", ("chromatin", "enriched"): "increased_accessibility_in",
    ("chromatin", "down"): "decreased_accessibility_in", ("chromatin", "NA"): "altered_accessibility_in",
    ("communication", "up"): "increased_communication_in", ("communication", "enriched"): "increased_communication_in",
    ("communication", "down"): "decreased_communication_in", ("communication", "exists"): "communicates_in",
    ("communication", "NA"): "altered_communication_in",
    ("spatial_colocalisation", "up"): "colocalises_in", ("spatial_colocalisation", "enriched"): "colocalises_in",
}
PRED_DIR = {"up": "increased_in", "down": "decreased_in", "enriched": "enriched_in", "depleted": "depleted_in",
            "exists": "exists_in", "associated_pos": "positively_associated_with",
            "associated_neg": "negatively_associated_with", "NA": "altered_in"}
NON_HUMAN = {"macaque": r"macaque|rhesus|primate", "mouse": r"\bmouse\b|\bmice\b|murine", "rat": r"\brats?\b|rodent"}
RATING = {}  # claims/store/claim_rating.tsv from 52_rating.py (filled in main when present)
HARMONISED_OF = {}  # vocab region label -> harmonised dataset regions it covers (filled in main)
PAPER_TO_DATASET = cfg.PAPER_TO_DATASET  # [claims] paper_to_dataset in config/project.toml


def split(v):
    if isinstance(v, list):
        return [str(x).strip() for x in v if str(x).strip()]
    return [x.strip() for x in re.split(r"[;,]", str(v or "")) if x.strip()]


VOCAB_COLS = ["label", "ontology", "term", "region_group", "harmonised_regions", "definition", "ontology_id",
              "ontology_label", "check"]
DATASET_COLS = ["dataset", "paper", "cohort", "title", "year", "doi", "accession", "access", "modality", "organism",
                "organism_id", "sample_type", "regions", "n_donors", "n_units", "donors_by_disease", "n_female", "n_male",
                "age_range", "pmi_range", "stage_available", "brain_banks", "country", "ancestry_summary",
                "ancestry_counts", "ancestry_basis", "ancestry_source", "disease_definitions", "staging", "sampling",
                "assay", "cohort_notes", "cohort_curated", "cohort_status"]


def vocab(name):
    f = VOCAB / f"{name}.tsv"
    if not f.exists():  # a project without vocabularies yet: every label is "not in vocabulary"
        return pd.DataFrame(columns=VOCAB_COLS).set_index("label")
    v = pd.read_csv(f, sep="\t").fillna("")
    for c in VOCAB_COLS:
        if c not in v:
            v[c] = ""
    return v.set_index("label")


def _read_tsv(path, columns, **kw):
    """A table that a new project may not have yet: empty with the expected columns."""
    if path.exists():
        return pd.read_csv(path, sep="\t", **kw)
    return pd.DataFrame(columns=columns)


def git_commit():
    try:
        return subprocess.run(["git", "-C", str(ROOT), "rev-parse", "--short", "HEAD"], capture_output=True,
                              text=True).stdout.strip()
    except OSError:
        return ""


# ---------------------------------------------------------------- datasets and cohorts
def load_datasets(papers):
    # inventory/cohorts.yaml = reviewed; inventory/cohorts_extracted.yaml = agent draft, used until reviewed
    reviewed, draft = ROOT / "inventory" / "cohorts.yaml", ROOT / "inventory" / "cohorts_extracted.yaml"
    cohorts_file = reviewed if reviewed.exists() else draft
    cohort_status = "reviewed" if reviewed.exists() else "draft (unreviewed extraction)"
    cohorts = {c["dataset"]: c for c in ((yaml.safe_load(cohorts_file.read_text()) or [])
                                          if cohorts_file.exists() else [])}
    rows, groups, units = [], [], []
    for f in sorted(H.glob("*_donors.tsv")):
        ds = f.name.replace("_donors.tsv", "")
        d = pd.read_csv(f, sep="\t", low_memory=False)
        donors = d.drop_duplicates("donor")
        p = papers.loc[ds] if ds in papers.index else pd.Series(dtype=object)
        c = cohorts.get(ds, {})
        anc = c.get("ancestry") or {}
        by_dis = donors.disease.astype(str).value_counts()
        rows.append(dict(
            dataset=ds, paper=ds, cohort=str(p.get("cohort", "") or ds), title=str(p.get("title", "")),
            year=p.get("year", ""), doi=str(p.get("doi", "")), accession=str(p.get("accession", "")),
            access=str(p.get("data_access", "")), modality=str(p.get("modality", "")),
            organism=c.get("organism", "Homo sapiens"), organism_id="NCBITaxon:9606",
            sample_type=c.get("sample_type", SAMPLE),
            regions=";".join(sorted(d.region.astype(str).unique())),
            n_donors=int(donors.donor.nunique()), n_units=len(d),
            donors_by_disease="; ".join(f"{k} {v}" for k, v in by_dis.items()),
            n_female=int((donors.sex == "F").sum()), n_male=int((donors.sex == "M").sum()),
            age_range=_rng(donors.age), pmi_range=_rng(donors.pmi),
            **{STAGE_COL: any(pd.to_numeric(donors[c], errors="coerce").notna().any()
                              for c in _contrasts.STAGE.get("columns", []) if c in donors)},
            brain_banks=_join(c.get("brain_banks")), country=str(c.get("country") or ""),
            ancestry_summary=str(anc.get("summary") or "not_reported"),
            ancestry_counts=json.dumps(anc.get("counts")) if anc.get("counts") else "",
            ancestry_basis=str(anc.get("basis") or "not_reported"), ancestry_source=str(anc.get("source") or ""),
            disease_definitions=_join(c.get("disease_definitions")),
            **{STAGING: _join(c.get(STAGING))}, sampling=_join(c.get("sampling")),
            assay=_join(c.get("assay")), cohort_notes=_join(c.get("notes")),
            cohort_curated=bool(c), cohort_status=cohort_status if c else "missing"))
        units.append(d[["dataset", "donor", "region", "disease"]].astype(str))  # in memory only, not stored
        for (reg, dis), g in d.groupby(["region", "disease"], observed=True):
            groups.append(dict(dataset=ds, region=reg, disease=dis, n_donors=g.donor.nunique()))
    units = pd.concat(units, ignore_index=True) if units else pd.DataFrame(columns=["dataset", "donor", "region", "disease"])
    empty_cols = [STAGE_COL if c == "stage_available" else STAGING if c == "staging" else c for c in DATASET_COLS]
    return (pd.DataFrame(rows, columns=None if rows else empty_cols),
            pd.DataFrame(groups, columns=None if groups else ["dataset", "region", "disease", "n_donors"]), units)


def _join(v):
    if isinstance(v, dict):
        return "; ".join(f"{k}: {_join(x)}" for k, x in v.items() if x)
    if isinstance(v, list):
        return "; ".join(_join(x) for x in v if x)
    return str(v or "")


def _rng(s):
    s = pd.to_numeric(s, errors="coerce").dropna()
    return f"{s.min():g}-{s.max():g}" if len(s) else ""


# ---------------------------------------------------------------- claim structure
def structure(claims, ctx, regions_v, cells_v, assays_v, problems):
    rows = []
    for _, c in claims.iterrows():
        x = ctx.get(c.id, {})
        direction = c.direction or "NA"
        null = str(x.get("null_claim", "")).lower() == "true"
        pred = "not_different_in" if null else PRED_TYPE.get((c.claim_type, direction), PRED_DIR.get(direction))
        # organism: human unless the claim's own evidence names another species
        blob = " ".join(str(c[k]) for k in ("statement", "comparison", "method", "reported_stat", "requires"))
        orgs = ["human"] + [o for o, pat in NON_HUMAN.items() if re.search(pat, blob, re.I)]
        # regions: as asserted; `multiple` resolved to the regions of the source data the claim was drawn from
        asserted = split(c.region)
        resolved = split(x.get("regions", "")) if (not asserted or "multiple" in asserted) else asserted
        bad = [r for r in resolved if r not in regions_v.index]
        if bad:
            problems.append((c.id, "region", f"not in vocabulary: {bad}"))
        cts = split(c.cell_type)
        badc = [t for t in cts if t not in cells_v.index]
        if badc:
            problems.append((c.id, "cell_type", f"not in vocabulary: {badc}"))
        if c.modality not in assays_v.index:
            problems.append((c.id, "modality", f"not in vocabulary: {c.modality}"))
        genes = split(c.genes)
        if c.cell_state:
            kind, label = "cell_state", c.cell_state
        elif c.claim_type in ("de_gene", "de_pathway", "chromatin") and genes:
            kind = "gene" if len(genes) == 1 else "gene_set"
            label = ", ".join(genes[:6]) + (f" (+{len(genes) - 6})" if len(genes) > 6 else "")
        elif c.claim_type == "genetic_enrichment":
            kind, label = "heritability", HERITABILITY
        else:
            kind, label = "cell_type", c.cell_type
        con = x.get("contrast", "")
        ckind, case, ref = CONTRASTS.get(con, ("unassigned", "", ""))
        if con not in CONTRASTS:
            problems.append((c.id, "contrast", f"code not in vocabulary: {con}"))
        rows.append(dict(
            claim_id=c.id, paper=c.paper, source_dataset=x.get("dataset", PAPER_TO_DATASET.get(c.paper, c.paper)),
            organisms=";".join(orgs), organism_ids=";".join(
                {"human": "NCBITaxon:9606", "macaque": "NCBITaxon:9544", "mouse": "NCBITaxon:10090",
                 "rat": "NCBITaxon:10116"}[o] for o in orgs),
            organism_basis=f"default (human {SAMPLE})" if len(orgs) == 1 else "keyword in claim evidence",
            organism_note="",
            sample_type=SAMPLE if orgs == ["human"] else f"{SAMPLE} + non-human data",
            region_asserted=";".join(asserted),
            regions=";".join(resolved), region_ids=";".join(regions_v.ontology_id.get(r, "") for r in resolved),
            region_groups=";".join(sorted({regions_v.region_group.get(r, "") for r in resolved} - {""})),
            cell_types=";".join(cts), cell_type_ids=";".join(cells_v.ontology_id.get(t, "") for t in cts),
            cell_state=c.cell_state, assay=c.modality, assay_id=assays_v.ontology_id.get(c.modality, ""),
            subject_kind=kind, subject=label, predicate=pred, direction=direction, null_claim=null,
            contrast=con, contrast_kind=ckind, case_group=case, reference_group=ref, comparison_text=c.comparison,
            claim_type=c.claim_type, centrality=c.centrality, unit_of_analysis=c.unit_of_analysis,
            test_template=c.test_template, testable_with=c.testable_with, context_status=x.get("status", ""),
            needs_review="; ".join(filter(None, [
                "non-human evidence" if len(orgs) > 1 else "",
                "region was 'multiple'" if "multiple" in asserted else "",
                "no contrast" if ckind in ("unassigned", "not_available") else ""]))))
    s = pd.DataFrame(rows)
    ov = ROOT / "claims" / "structure_overrides.tsv"
    if ov.exists():
        for _, o in pd.read_csv(ov, sep="\t", comment="#").fillna("").iterrows():
            m = s.claim_id.eq(o.claim_id)
            if not m.any() or o.field not in s:
                problems.append((o.claim_id, "override", f"unknown claim or field {o.field}"))
                continue
            s.loc[m, o.field] = o.value
            s.loc[m, "needs_review"] = s.loc[m, "needs_review"].str.replace(o.get("resolves", "") or "@@", "",
                                                                            regex=False).str.strip("; ")
    return s


# ---------------------------------------------------------------- evidence (claim x dataset)
def load_tests():
    parts = []
    for sub, role in (("internal", "source"), ("internal_gwas", "source"), ("cross", "replication"),
                      ("cross_gwas", "replication")):
        for f in sorted((RES / "tests" / sub).glob("*.tsv")):
            t = pd.read_csv(f, sep="\t", low_memory=False)
            t["role"], t["engine"] = role, "gwas" if sub.endswith("gwas") else "main"
            parts.append(t)
    keep = ["claim_id", "dataset", "role", "engine", "variant", "contrast", "region_match", "template_used",
            "outcome", "effect", "z", "p", "n_case", "n_ref", "regions_tested", "classes_tested", "set_test", "note"]
    t = pd.concat(parts, ignore_index=True) if parts else pd.DataFrame(columns=keep)  # no tests run yet
    for k in keep:
        if k not in t:
            t[k] = np.nan
    t = t[keep]
    t.loc[t.role.eq("source"), "region_match"] = "same"
    # one claim can have rows from both engines only if it is genetic; prefer the GWAS engine for those
    t = t.sort_values("engine").drop_duplicates(["claim_id", "dataset", "variant"], keep="last")
    return t


# ---------------------------------------------------------------- findings and questions (claimstore.store)
def load_findings(problems):
    """claims/findings.yaml -> claim rows, context rows and test rows (evidence resolved from results tables)."""
    from claimstore.store import FINDINGS, finding_tests, resolve_evidence
    items = yaml.safe_load(FINDINGS.read_text()) if FINDINGS.exists() else []
    rows, ctx, tests = [], {}, []
    for f in items or []:
        ev = resolve_evidence(f["evidence"])
        if ev.empty:
            problems.append((f["id"], "evidence", "no rows matched the evidence reference"))
        regions = sorted({r for v in ev.get("regions", pd.Series(dtype=str)) for r in str(v).split(";") if r and r != "nan"})
        if not regions and len(ev):  # results without a region column: the regions the datasets sampled
            regions = sorted({r for d in ev.dataset if (H / f"{d}_donors.tsv").exists()
                              for r in pd.read_csv(H / f"{d}_donors.tsv", sep="\t", usecols=["region"]).region.astype(str)} - {"other"})
        region = f["region"] if f["region"] != "multiple" else "multiple"
        rows.append({**{k: f.get(k, "") for k in ("id", "paper", "statement", "quote", "location", "centrality", "modality",
                                                  "cell_type", "cell_state", "claim_type", "direction", "comparison", "genes",
                                                  "gene_source", "reported_stat", "unit_of_analysis", "method", "testable_with",
                                                  "test_template", "requires", "notes")},
                     "region": region, "origin": "finding", "themes_f": f.get("themes", ""),
                     "analysis_script": f["analysis"].get("script", ""), "analysis_results": f["analysis"].get("results", ""),
                     "analysis_date": f["analysis"].get("date", "")})
        ctx[f["id"]] = {"dataset": "", "regions": ";".join(regions) or f["region"], "contrast": f["contrast"],
                        "null_claim": str(bool(f.get("null_claim"))), "status": "ok"}
        tests += finding_tests(f, ev)
    return pd.DataFrame(rows), ctx, pd.DataFrame(tests)


def load_background():
    bf = ROOT / "claims" / "background.yaml"
    items = (yaml.safe_load(bf.read_text()) or []) if bf.exists() else []
    cols = ["id", "statement", "evidence_type", "references", "statements", "establishes", "supports_claims", "consistent_with", "contrasts_with", "regions",
            "measurement_note", "added"]
    return pd.DataFrame([{k: (json.dumps(b.get(k)) if isinstance(b.get(k), (list, dict)) else b.get(k, "")) for k in cols}
                         for b in items], columns=cols)


def load_questions():
    qf = ROOT / "claims" / "questions.yaml"
    items = (yaml.safe_load(qf.read_text()) or []) if qf.exists() else []
    cols = ["id", "text", "asked_by", "date", "status", "answer", "records", "proposed"]
    return pd.DataFrame([{k: (json.dumps(q.get(k)) if isinstance(q.get(k), (list, dict)) else q.get(k, "")) for k in cols}
                         for q in items], columns=cols)


# ---------------------------------------------------------------- relations
def relations(struct, claims):
    rel = []
    cs = _read_tsv(ROOT / "claims" / "consensus_statements.tsv", ["id", "statement", "short", "our_evidence", "verdict",
                                                                    "supporting_claims", "opposing_claims"]).fillna("")
    stmts = cs[["id", "statement", "short", "our_evidence", "verdict"]].rename(columns={"id": "statement_id"})
    for _, r in cs.iterrows():
        for col, stance in (("supporting_claims", "supports"), ("opposing_claims", "opposes")):
            for cid in split(r[col]):
                rel.append(dict(claim_a=cid, relation=f"{stance}_statement", claim_b=r.id, basis="curated",
                                note=r.short))
    # computed: composition claims about the same cell type/state, same contrast family, different papers
    comp = struct[struct.claim_type.eq("composition") & struct.contrast_kind.isin(["case_control", "stage"])]
    sign = comp.predicate.map({"increased_proportion_in": 1, "decreased_proportion_in": -1})
    comp = comp.assign(sign=sign, key=comp.cell_types + "|" + comp.cell_state.str.lower())
    for _, g in comp.dropna(subset=["sign"]).groupby("key"):
        g = g.reset_index(drop=True)
        for i in range(len(g)):
            for j in range(i + 1, len(g)):
                a, b = g.loc[i], g.loc[j]
                if a.paper == b.paper:
                    continue
                rel.append(dict(claim_a=a.claim_id, relation="agrees_with" if a.sign == b.sign else "conflicts_with",
                                claim_b=b.claim_id, basis="computed",
                                note=f"composition of {a.cell_types} {a.cell_state}".strip()))
    return pd.DataFrame(rel, columns=["claim_a", "relation", "claim_b", "basis", "note"]), stmts


# ---------------------------------------------------------------- contextualised text
def _i(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return 0


def context_text(s, ds, summ, ev, units):
    d = ds.get(s.source_dataset, {})
    org = ("Human" if s.organisms == "human"
           else f"Human (+ {s.organisms.replace('human;', '').replace(';', ', ')} data)")
    regs = s.regions.split(";") if s.regions else []
    region = (f"multiple regions (re-tested in: {', '.join(regs)})" if "multiple" in s.region_asserted
              else ", ".join(regs) or "region not specified")
    bits = [f"{org} {SAMPLE}", region,
            s.cell_types.replace(";", ", ") + (f" / {s.cell_state}" if s.cell_state else ""), s.assay]
    if s.get("organism_note"):
        bits.append(s.organism_note)
    if d.get("sampling"):
        bits.append(d["sampling"])
    groups = f"{s.case_group} vs {s.reference_group}" if s.case_group else s.comparison_text
    # donors of the source dataset in the claim's regions (all regions if none match)
    g = units[units.dataset.eq(s.source_dataset)]
    want = set(regs) | {h for r in regs for h in split(HARMONISED_OF.get(r, ""))}  # e.g. a region group -> its sub-regions
    g = g[g.region.isin(want)] if g.region.isin(want).any() else g
    by = g.drop_duplicates("donor").disease.value_counts()
    donors = "; ".join(f"{k} {v}" for k, v in by.items())
    cohort = (f"cohort {d.get('cohort', s.source_dataset)} (donors in these regions: {donors}); "
              f"ancestry: {d.get('ancestry_summary', 'not_reported')}")
    if d.get("brain_banks"):
        cohort += f"; brain bank: {d['brain_banks']}"
    if not s.source_dataset:  # a finding of this study: evidence is the donor-level analysis across datasets
        own = ev[ev.claim_id.eq(s.claim_id)]
        cohort = (f"this study's donor-level analysis of {own.dataset.nunique()} datasets "
                  f"({int(own.n_case.fillna(0).sum())} case and {int(own.n_ref.fillna(0).sum())} reference donor units)")
    sm = summ.get(s.claim_id, {})
    src = ev[ev.claim_id.eq(s.claim_id) & ev.role.eq("source")].outcome.dropna()
    evid = (f"Re-test in source data: {sm.get('own_data_outcome') or (src.iloc[0] if len(src) else 'not run')}. "
            f"Other datasets: tested in {_i(sm.get('k'))} ({_i(sm.get('n_reproduced'))} reproduced, "
            f"{_i(sm.get('n_direction_only'))} direction only, {_i(sm.get('n_not_reproduced'))} not reproduced, "
            f"{_i(sm.get('n_contradicted'))} contradicted; same region: {_i(sm.get('same_region_rep'))} of "
            f"{_i(sm.get('same_region_k'))} reproduced).")
    rt = RATING.get(s.claim_id)
    if not s.source_dataset:
        evid = ""
    if rt:
        eff = (f" Pooled effect (approx. SMD, positive = as claimed) {rt['est']:.2f}, 95% CI {rt['lo']:.2f} to {rt['hi']:.2f}."
               if pd.notna(rt.get("est")) else "")
        evid += (f" Rating: {rt['verdict']}, {rt['certainty'].lower()} certainty ({rt['k_cohorts']} independent cohorts)."
                 f"{eff} Concerns: {rt['reasons']}.")
    return ("Context: " + "; ".join(bits) + f". Comparison: {groups} (unit in source analysis: "
            f"{s.unit_of_analysis}). Source: {cohort}."), evid


def main():
    OUT_GIT.mkdir(parents=True, exist_ok=True)
    OUT_DB.mkdir(parents=True, exist_ok=True)
    problems = []
    claims = pd.DataFrame(yaml.safe_load((ROOT / "claims" / "claims.yaml").read_text())).fillna("")
    claims["origin"] = "published"
    ctx = pd.read_csv(ROOT / "claims" / "contexts.tsv", sep="\t", dtype=str).fillna("").set_index("claim_id")
    ctx = ctx.to_dict("index")
    fclaims, fctx, ftests = load_findings(problems)
    claims = pd.concat([claims, fclaims], ignore_index=True).fillna("")
    ctx.update(fctx)
    papers = pd.read_csv(ROOT / "inventory" / "papers.tsv", sep="\t").fillna("").set_index("paper_id")
    regions_v, cells_v, assays_v = vocab("regions"), vocab("cell_types"), vocab("assays")
    rating_f = OUT_GIT / "claim_rating.tsv"
    rating = pd.read_csv(rating_f, sep="\t") if rating_f.exists() else pd.DataFrame(columns=["claim_id"])
    RATING.update(rating.set_index("claim_id").to_dict("index"))
    HARMONISED_OF.update(regions_v.harmonised_regions.to_dict())

    datasets, donor_groups, donor_units = load_datasets(papers)
    for _, d in datasets.iterrows():
        if not d.cohort_curated:
            problems.append((d.dataset, "cohort", "no cohort entry (inventory/cohorts.yaml or cohorts_extracted.yaml)"))
        elif d.ancestry_basis == "not_reported":
            problems.append((d.dataset, "ancestry", "not reported"))
    struct = structure(claims, ctx, regions_v, cells_v, assays_v, problems)
    is_finding = struct.claim_id.isin(fclaims.id) if len(fclaims) else pd.Series(False, index=struct.index)
    struct.loc[is_finding, "needs_review"] = struct.loc[is_finding, "needs_review"].str.replace("region was 'multiple'", "").str.strip("; ")
    struct.loc[is_finding, "source_dataset"] = ""
    for cid in struct.claim_id[~struct.source_dataset.isin(datasets.dataset) & ~is_finding]:
        problems.append((cid, "source_dataset", "source dataset not harmonised"))
    tests = load_tests()
    tests["d_native"], tests["se_native"] = np.nan, np.nan
    if len(ftests):
        tests = pd.concat([tests, ftests], ignore_index=True)
    rel, stmts = relations(struct, claims)
    ev_cols = ["claim_id", "rank", "tier", "kind", "themes", "own_data_outcome", "k", "n_reproduced",
               "n_direction_only", "n_not_reproduced", "n_contradicted", "share_reproducing", "z", "p", "I2",
               "same_region_k", "same_region_rep", "other_region_k", "other_region_rep", "robustness"]
    ev_sum = _read_tsv(ROOT / "claims" / "claim_evidence.tsv", ev_cols, low_memory=False)
    ev_sum = ev_sum[["claim_id", "rank", "tier", "kind", "themes", "own_data_outcome", "k", "n_reproduced",
                     "n_direction_only", "n_not_reproduced", "n_contradicted", "share_reproducing", "z", "p", "I2",
                     "same_region_k", "same_region_rep", "other_region_k", "other_region_rep", "robustness"]]
    sigs = _read_tsv(ROOT / "claims" / "signatures_all.tsv.gz", ["signature_id", "paper", "claim_ids", "kind", "cell_type",
                     "cell_state", "comparison", "gene", "rank", "stat", "padj", "source", "selection", "variant",
                     "gene_reported"], low_memory=False)
    for c in ["comparison", "rank", "stat", "padj", "selection"]:  # optional in a project's own signature table
        if c not in sigs:
            sigs[c] = np.nan
    sig_meta = sigs.drop_duplicates("signature_id")[["signature_id", "paper", "kind", "cell_type", "cell_state",
                                                    "comparison", "source", "selection"]]
    sig_claim = (sigs.drop_duplicates("signature_id")[["signature_id", "claim_ids"]]
                 .assign(claim_id=lambda x: x.claim_ids.str.split(";")).explode("claim_id")[["claim_id", "signature_id"]].dropna())
    sig_gene = sigs[["signature_id", "gene", "rank", "stat", "padj"]]
    disc = _read_tsv(ROOT / "claims" / "source_discrepancies.tsv", ["claim_id", "finding", "evidence", "found_by"]).fillna("")
    sources = papers.reset_index()[["paper_id", "first_author", "year", "title", "journal_or_preprint", "doi", "pmid",
                                    "modality", "regions", "screening", "cohort", "reuses_data_of"]]
    sources.insert(1, "kind", "paper")

    # claims table: original fields + structure + contextualised text
    summ = ev_sum.set_index("claim_id").to_dict("index")
    dsd = datasets.set_index("dataset").to_dict("index")
    ctx_txt = struct.apply(lambda s: context_text(s, dsd, summ, tests, donor_units), axis=1, result_type="expand")
    struct["context_text"], struct["evidence_text"] = ctx_txt[0], ctx_txt[1]
    for col in ("origin", "themes_f", "analysis_script", "analysis_results", "analysis_date"):
        if col not in claims:
            claims[col] = ""
    base = claims[["id", "statement", "quote", "location", "comparison", "genes", "gene_source", "reported_stat",
                   "method", "requires", "notes", "origin", "themes_f", "analysis_script", "analysis_results",
                   "analysis_date"]].rename(columns={"id": "claim_id", "comparison": "comparison_raw"})
    base["genes"] = base.genes.map(lambda g: ";".join(split(g)))
    full = struct.merge(base, on="claim_id")
    full["source_kind"] = np.where(full.origin.eq("finding"), "our_analysis", "paper")
    full["source_id"] = np.where(full.origin.eq("finding"), full.analysis_script, full.paper)
    full["status"] = np.where(full.origin.eq("finding"), "rated", "asserted")
    full = full.merge(ev_sum[["claim_id", "tier", "rank", "themes"]], on="claim_id", how="left")
    full["themes"] = full.themes.where(full.themes.notna() & full.themes.ne(""), full.themes_f)
    full = full.drop(columns=["themes_f"])
    rcols = ["verdict", "certainty", "rating", "family", "start", "k_cohorts", "est", "lo", "hi", "i2", "concordance",
             "z_comb", "detect", "donors", "rob", "inconsistency", "indirectness", "imprecision", "dependence", "upgrade",
             "reasons", "not_assessed", "scope_tested", "scope_held", "caveats"]
    full = full.merge(rating[["claim_id"] + [c for c in rcols if c in rating]], on="claim_id", how="left")
    for c in rcols:  # not rated yet: the columns exist, empty
        if c not in full:
            full[c] = None

    # ---- write git outputs
    struct.drop(columns=["context_text", "evidence_text"]).to_csv(OUT_GIT / "claim_structure.tsv", sep="\t",
                                                                   index=False)
    with open(OUT_GIT / "claims_contextualised.jsonl", "w") as fh:
        for _, r in full.iterrows():
            fh.write(json.dumps(dict(
                id=r.claim_id, source=r.paper, statement=r.statement, context=r.context_text,
                evidence=r.evidence_text, subject=r.subject, predicate=r.predicate, contrast=r.contrast,
                organisms=r.organisms, regions=r.regions, region_ids=r.region_ids, cell_types=r.cell_types,
                cell_type_ids=r.cell_type_ids, assay=r.assay, source_dataset=r.source_dataset,
                verdict=r.get("verdict"), certainty=r.get("certainty"), origin=r.origin, location=r.location, quote=r.quote),
                ensure_ascii=False) + "\n")
    pd.DataFrame(problems, columns=["item", "field", "problem"]).to_csv(OUT_GIT / "validation.tsv", sep="\t",
                                                                         index=False)

    # ---- sqlite
    db = OUT_DB / "claim_store.sqlite"
    db.unlink(missing_ok=True)
    con = sqlite3.connect(db)
    tables = dict(claim=full, dataset=datasets, donor_group=donor_groups, source=sources,
                  claim_dataset=tests, claim_summary=ev_sum, claim_relation=rel, consensus_statement=stmts,
                  signature=sig_meta, claim_signature=sig_claim, signature_gene=sig_gene,
                  source_discrepancy=disc, question=load_questions(), background=load_background())
    for k in ("regions", "cell_types", "diseases", "assays", "organisms"):
        tables[f"vocab_{k}"] = _read_tsv(VOCAB / f"{k}.tsv", VOCAB_COLS)
    tables["vocab_contrasts"] = pd.DataFrame([dict(code=k, kind=v[0], case_group=v[1], reference_group=v[2])
                                              for k, v in CONTRASTS.items()],
                                             columns=["code", "kind", "case_group", "reference_group"])
    tables["build_info"] = pd.DataFrame([dict(commit=git_commit(), built=pd.Timestamp.now().isoformat(),
                                              n_claims=len(full), n_datasets=len(datasets))])
    for name, df in tables.items():
        df.to_sql(name, con, index=False)
    for idx in ("claim_dataset(claim_id)", "claim_dataset(dataset)", "claim_signature(claim_id)",
                "signature_gene(gene)", "claim_relation(claim_a)", "claim_relation(claim_b)"):
        con.execute(f"CREATE INDEX IF NOT EXISTS ix_{idx.replace('(', '_').rstrip(')')} ON {idx}")
    con.commit()
    con.close()
    print(f"claims {len(full)}, datasets {len(datasets)}, claim_dataset {len(tests)}, relations {len(rel)}, "
          f"problems {len(problems)} -> {db}")
    print(pd.DataFrame(problems, columns=["item", "field", "problem"]).field.value_counts().to_string())


if __name__ == "__main__":
    main()

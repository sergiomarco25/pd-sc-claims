"""Evidence rating v2: verdict + certainty (GRADE-style) + scope of validity, for every claim in the store.

Rules: claims/RATING.md (version RATING_VERSION below; keep the two equal; tests/test_rating_doc.py checks it).
Facts about the project's datasets and test templates come from the [claims] and [rating] tables of its config file
(claimstore.config): which templates form each family, datasets of sorted nuclei, datasets that may share donors, papers
whose claims were made on another cohort, and what the rating does not yet assess.

Inputs: results/store/claim_store.sqlite (51b: claim, claim_dataset, dataset), claims/signature_coverage.tsv,
        claims/claim_evidence.tsv (old tiers, for comparison only).
Outputs: claims/store/claim_rating.tsv        one row per claim: verdict, certainty, every domain with its reason
         claims/store/rating_vs_tiers.tsv     cross-tabulation old tier x new verdict/certainty
         claims/store/rating_tests.tsv        per test: oriented z, approximate SMD, SE, weight, flags

Evidence families
  effect   group contrasts and correlations (composition, DE genes, gene-set and state-score shifts, trajectories).
           Each test's z is oriented (positive = in the claimed direction, from the engine's outcome) and converted
           to an approximate standardised mean difference: d = z * sqrt(1/n1 + 1/n2) / sqrt(m) for two groups
           (m = regions x classes combined by the engine's Stouffer step), or via r = z / sqrt(n) for correlations.
           Datasets pooled by DerSimonian-Laird random effects.
  pattern  claims that a structure exists or that two contrasts agree (state existence, concordance between two
           diseases, genetic enrichment). Their z comes from permutation or enrichment tests, not donor contrasts, so the
           verdict uses detection across datasets (weighted Stouffer of oriented z and the share detecting it).
Thresholds are initial values (ROADMAP 3.1-3.2) and are collected in THRESH for the sensitivity analysis.

Run through the project's entry point, e.g. python -I scripts/52_rating.py [--db claim_store.sqlite] [--out folder]
     (defaults: the configured store and claims/store/; the release builder re-rates a public store with --db/--out)
"""
import re
import sqlite3
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from claimstore import config as cfg

ROOT = cfg.REPO

DB = cfg.STORE_DB
OUT = ROOT / "claims" / "store"
RATING_VERSION = "2.0.0"  # claims/RATING.md; written to every row of claim_rating.tsv
THRESH = dict(equiv=0.2, trivial=0.2, important=0.5, large=0.8, i2_serious=0.5, i2_very=0.75, conc_serious=0.7,
              conc_very=0.6, donors_ois=200, donors_min=30, weight_major=0.5, detect_supported=2 / 3)
LEVELS = ["Very low", "Low", "Moderate", "High"]
# project facts ([rating] and [claims] of the config file)
EFFECT_T = set(cfg.RATING.get("effect_templates", []))
PATTERN_T = set(cfg.RATING.get("pattern_templates", []))
ENRICHED_SORT = set(cfg.RATING.get("enriched_sort", []))  # sorted nuclei: composition is not tissue composition
# datasets that may share donors (not verified); reported as a caveat, not merged
SHARED_BANK = {k: set(v) for k, v in cfg.RATING.get("shared_banks", {}).items()}
SHARED_CAVEAT = cfg.RATING.get("shared_bank_caveat", "{dataset} may share donors with {others} (not verified)")
NOT_ASSESSED = cfg.RATING.get("not_assessed", "")
PAPER_TO_DATASET = cfg.PAPER_TO_DATASET


def dl_meta(y, se):
    """DerSimonian-Laird random effects: pooled estimate, SE, tau2, I2, 95% prediction interval."""
    y, v = np.asarray(y, float), np.asarray(se, float) ** 2
    w = 1 / v
    fe = (w * y).sum() / w.sum()
    q = (w * (y - fe) ** 2).sum()
    k = len(y)
    c = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 and c > 0 else 0.0
    ws = 1 / (v + tau2)
    est = (ws * y).sum() / ws.sum()
    se_est = np.sqrt(1 / ws.sum())
    i2 = max(0.0, (q - (k - 1)) / q) if k > 1 and q > 0 else 0.0
    if k >= 3:
        t = stats.t.ppf(0.975, k - 2)
        pi = (est - t * np.sqrt(tau2 + se_est ** 2), est + t * np.sqrt(tau2 + se_est ** 2))
    else:
        pi = (np.nan, np.nan)
    return dict(est=est, se=se_est, tau2=tau2, i2=i2, pi_lo=pi[0], pi_hi=pi[1], weights=ws / ws.sum())


def orient(row, expected):
    """Oriented z: positive = in the claimed direction. The engine's outcome already encodes the orientation."""
    z = abs(float(row.z))
    if expected == 0:
        return float(row.z)  # direction-free claims keep the data sign
    return z if row.outcome in ("reproduced", "direction_only") else -z


def to_smd(zo, n1, n2, m):
    m = max(1, int(m) if np.isfinite(m) else 1)
    if n2 and n2 > 0 and n1 and n1 > 0:
        d = zo * np.sqrt(1 / n1 + 1 / n2) / np.sqrt(m)
        d = float(np.clip(d, -3, 3))
        se = np.sqrt(1 / n1 + 1 / n2 + d ** 2 / (2 * (n1 + n2)))
        return d, se, "two-group"
    n = n1 if n1 else np.nan
    if not n or n < 5:
        return np.nan, np.nan, "too few units"
    r = float(np.clip(zo / np.sqrt(n) / np.sqrt(m), -0.95, 0.95))
    d = 2 * r / np.sqrt(1 - r ** 2)
    se = np.sqrt(4 / ((n - 1) * (1 - r ** 2)))
    return d, se, "correlation"


def load():
    con = sqlite3.connect(DB)
    claims = pd.read_sql("SELECT * FROM claim", con)
    tests = pd.read_sql("SELECT * FROM claim_dataset", con)
    ds = pd.read_sql("SELECT dataset, cohort FROM dataset", con)
    con.close()
    for c in ("z", "p", "effect", "n_case", "n_ref", "d_native", "se_native"):
        if c not in tests:
            tests[c] = np.nan
        tests[c] = pd.to_numeric(tests[c], errors="coerce")
    raw = []  # n_tests is not in the store: read it from the engine outputs
    base = cfg.RESULTS / "tests"
    for sub in ("internal", "internal_gwas", "cross", "cross_gwas"):
        for f in sorted((base / sub).glob("*.tsv")):
            t = pd.read_csv(f, sep="\t", usecols=lambda c: c in ("claim_id", "dataset", "variant", "n_tests"))
            raw.append(t)
    nt = (pd.concat(raw).drop_duplicates(["claim_id", "dataset", "variant"]) if raw  # none yet in a new project
          else pd.DataFrame(columns=["claim_id", "dataset", "variant", "n_tests"]))
    nt["variant"] = nt["variant"].astype(object)  # all-empty in a project without claim variants: read as float
    tests = tests.merge(nt, on=["claim_id", "dataset", "variant"], how="left")
    cov_f, old_f = ROOT / "claims" / "signature_coverage.tsv", ROOT / "claims" / "claim_evidence.tsv"
    cov = (pd.read_csv(cov_f, sep="\t").set_index("claim_id").coverage if cov_f.exists()
           else pd.Series(dtype=object, name="coverage"))
    old = (pd.read_csv(old_f, sep="\t", low_memory=False).set_index("claim_id").tier if old_f.exists()
           else pd.Series(dtype=object, name="tier"))
    return claims, tests, ds, cov, old


def pick_variant(t):
    order = {"both": 0, "revised": 1, "preprint": 2}
    t = t.assign(_o=t.variant.map(order).fillna(3))
    return t.sort_values("_o").drop_duplicates(["claim_id", "dataset"]).drop(columns="_o")


def rate_claim(c, t, cohort_of, cov):
    expected = {"up": 1, "enriched": 1, "associated_pos": 1, "exists": 1, "down": -1, "depleted": -1,
                "associated_neg": -1}.get(c.direction, 0)
    source = PAPER_TO_DATASET.get(c.paper, c.source_dataset)
    inf = t[t.outcome.isin(["reproduced", "direction_only", "not_reproduced", "contradicted"]) & t.z.notna()].copy()
    out = dict(claim_id=c.claim_id, family="", verdict="Untested", certainty="Very low", start="Very low",
               k_tests=len(inf), k_cohorts=0, est=np.nan, lo=np.nan, hi=np.nan, i2=np.nan, concordance=np.nan,
               z_comb=np.nan, detect=np.nan,
               donors=0, rob=0, inconsistency=0, indirectness=0, imprecision=0, dependence=0, upgrade=0,
               reasons="", not_assessed="", scope_tested="", scope_held="", caveats="")
    if inf.empty:
        out["reasons"] = "no informative test; rests on the paper's own analysis" + (
            " at cell level" if "cell" in str(c.unit_of_analysis).lower() else "")
        return out, []
    inf["zo"] = inf.apply(lambda r: orient(r, expected), axis=1)
    templates = set(inf.template_used)
    family = "pattern" if templates & PATTERN_T and not templates & EFFECT_T else "effect"
    out["family"] = family
    inf["cohort"] = inf.dataset.map(cohort_of)
    out["k_cohorts"] = inf.cohort.nunique()
    out["donors"] = int((inf.n_case.fillna(0) + inf.n_ref.fillna(0)).sum())
    reasons, flags_rows = [], []

    # ---- per-test flags
    claimed_regions = set(str(c.regions).split(";")) if c.regions else set()
    for i, r in inf.iterrows():
        f = []
        if r.role != "source" and r.region_match == "other":
            f.append("other region")
        if c.cell_state and str(cov.get(c.claim_id, "")) in ("partial", "text_only", "none") \
                and r.template_used in ("state_genes_shift", "geneset_shift", "state_existence"):
            f.append("partial proxy for the state")
        if r.template_used == "composition_clr" and r.dataset in ENRICHED_SORT:
            f.append("sorted nuclei")
        inf.at[i, "flag"] = ";".join(f)

    # ---- verdict and pooled estimate
    if family == "effect":
        conv = inf.apply(lambda r: (float(r.d_native) * (1 if expected == 0 or r.zo * np.sign(r.d_native) >= 0 else -1),
                                    float(r.se_native), "native") if pd.notna(r.get("d_native")) else
                         to_smd(r.zo, r.n_case, r.n_ref, r.n_tests if pd.notna(r.n_tests) else 1), axis=1)
        inf["d"], inf["se_d"], inf["conv"] = [x[0] for x in conv], [x[1] for x in conv], [x[2] for x in conv]
        ok = inf[np.isfinite(inf.d) & np.isfinite(inf.se_d)]
        if ok.empty:
            out["reasons"] = "tests too small to estimate an effect"
            return out, []
        m = dl_meta(ok.d, ok.se_d)
        inf.loc[ok.index, "weight"] = m["weights"]
        est, se = m["est"], m["se"]
        lo, hi = est - 1.96 * se, est + 1.96 * se
        out.update(est=est, lo=lo, hi=hi, i2=m["i2"])
        if c.null_claim in (1, True, "True", "1"):
            out["verdict"] = "Supported (no effect)" if -THRESH["equiv"] < lo and hi < THRESH["equiv"] else \
                "Contradicted" if (lo > 0 or hi < 0) else "Inconclusive"
        elif expected == 0:
            out["verdict"] = "Supported" if (lo > 0 or hi < 0) else ("No effect" if -THRESH["equiv"] < lo and hi < THRESH["equiv"] else "Inconclusive")
        else:
            out["verdict"] = "Supported" if lo > 0 else "Contradicted" if hi < 0 else \
                "No effect" if -THRESH["equiv"] < lo and hi < THRESH["equiv"] else "Inconclusive"
        weights = inf.weight.fillna(0)
        pi_cross = m["pi_lo"] < 0 < m["pi_hi"] if np.isfinite(m["pi_lo"]) else False
    else:
        w = np.sqrt(inf.n_case.fillna(1).clip(lower=1) + inf.n_ref.fillna(0))
        zc = float((w * inf.zo).sum() / np.sqrt((w ** 2).sum()))
        detect = float(inf.outcome.eq("reproduced").mean())
        same_dir = float((inf.zo > 0).mean())
        out.update(z_comb=zc, detect=detect, i2=np.nan)
        if c.null_claim in (1, True, "True", "1"):
            out["verdict"] = "Supported (no effect)" if detect >= THRESH["detect_supported"] else "Contradicted" \
                if inf.outcome.eq("contradicted").mean() > .5 else "Inconclusive"
        elif detect >= THRESH["detect_supported"] or (zc > 1.96 and detect >= .5):
            out["verdict"] = "Supported"
        elif zc < -1.96 or inf.outcome.eq("contradicted").mean() > .5:
            out["verdict"] = "Contradicted"
        elif same_dir >= THRESH["detect_supported"]:
            out["verdict"] = "Partly supported"  # present in the claimed direction, but not meeting the claim's
            #                                      specific criterion (e.g. a programme, not a distinct subpopulation)
        else:
            out["verdict"] = "Inconclusive"
        weights = w / w.sum()
        inf["weight"] = weights
        pi_cross = False
    conc = float((inf.zo > 0).mean()) if expected != 0 else np.nan  # share in the claimed direction (reported)
    out["concordance"] = conc
    # agreement among datasets is judged against the pooled direction, not the claimed one: a contradicted claim
    # whose datasets agree on the opposite direction is consistent evidence
    pooled_sign = np.sign(out["est"]) if family == "effect" and np.isfinite(out["est"]) else 1
    agree = float((np.sign(inf.zo) == pooled_sign).mean()) if expected != 0 else np.nan

    # ---- starting level
    k = out["k_cohorts"]
    out["start"] = "High" if k >= 3 else "Moderate" if k == 2 else "Low"
    level = LEVELS.index(out["start"])

    # ---- risk of bias
    w_bias = float(weights[inf["flag"].str.contains("sorted nuclei", na=False)].sum())
    if w_bias > THRESH["weight_major"]:
        out["rob"] = -1
        reasons.append(f"risk of bias: {w_bias:.0%} of weight from sorted nuclei")
    # ---- inconsistency
    if family == "effect" and len(inf) >= 2:
        i2 = out["i2"]
        if i2 > THRESH["i2_very"] and np.isfinite(agree) and agree < THRESH["conc_very"]:
            out["inconsistency"] = -2
            reasons.append(f"inconsistency: I² {i2:.0%}, {agree:.0%} of datasets in the pooled direction")
        elif i2 >= THRESH["i2_serious"] or (np.isfinite(agree) and agree < THRESH["conc_serious"]) or \
                (pi_cross and out["verdict"] in ("Supported", "Contradicted")):
            out["inconsistency"] = -1
            reasons.append(f"inconsistency: I² {i2:.0%}" + (f", {agree:.0%} of datasets in the pooled direction" if np.isfinite(agree) else "")
                           + (", prediction interval crosses 0" if pi_cross else ""))
    elif family == "pattern" and len(inf) >= 2:
        detect = float(inf.outcome.eq("reproduced").mean())
        if detect < .5 and out["verdict"].startswith(("Supported", "Partly")):
            out["inconsistency"] = -1
            reasons.append(f"inconsistency: detected in {detect:.0%} of datasets")
    # ---- indirectness
    ind = inf["flag"].fillna("")
    w_ind = float(weights[ind.str.len() > 0].sum())
    w_two = float(weights[ind.str.count(";") >= 1].sum())
    nonhuman = c.organisms != "human" and "macaque" in str(c.organisms)
    if nonhuman or w_two > THRESH["weight_major"]:
        out["indirectness"] = -2
        reasons.append("indirectness: " + ("non-human evidence in the claim" if nonhuman else
                                           f"{w_two:.0%} of weight with two mismatches (region, state proxy)"))
    elif w_ind > THRESH["weight_major"]:
        out["indirectness"] = -1
        kinds = sorted({k for f in ind for k in f.split(";") if k})
        reasons.append(f"indirectness: {w_ind:.0%} of weight with " + ", ".join(kinds))
    # ---- imprecision
    donors = out["donors"]
    if len(inf) == 1 or donors < THRESH["donors_min"]:
        out["imprecision"] = -2
        reasons.append(f"imprecision: {len(inf)} test(s), {donors} donor units")
    elif family == "effect":
        lo_a, hi_a = sorted([abs(out["lo"]), abs(out["hi"])]) if (out["lo"] > 0 or out["hi"] < 0) else (0, max(abs(out["lo"]), abs(out["hi"])))
        wide = lo_a < THRESH["trivial"] and hi_a > THRESH["important"]
        if wide or donors < THRESH["donors_ois"]:
            out["imprecision"] = -1
            reasons.append("imprecision: " + (f"CI {out['lo']:.2f} to {out['hi']:.2f} spans trivial and important effects"
                                              if wide else f"{donors} donor units (< {THRESH['donors_ois']})"))
    elif donors < THRESH["donors_ois"]:
        out["imprecision"] = -1
        reasons.append(f"imprecision: {donors} donor units")
    # ---- dependence
    src_cohort = cohort_of.get(source, source)
    w_src = float(weights[inf.cohort.eq(src_cohort)].sum())
    if w_src >= 0.999:
        out["dependence"] = -2
        reasons.append("dependence: all evidence from the source cohort")
    elif w_src > THRESH["weight_major"]:
        out["dependence"] = -1
        reasons.append(f"dependence: {w_src:.0%} of weight from the source cohort")
    caveats = []
    for anchor in sorted(SHARED_BANK):
        shared = sorted(set(inf.dataset) & SHARED_BANK[anchor]) if anchor in set(inf.dataset) else []
        if shared:
            caveats.append(SHARED_CAVEAT.format(dataset=anchor, others=", ".join(shared)))
    if caveats:
        out["caveats"] = "; ".join(caveats)
    # ---- upgrades
    worst = min(out["rob"], out["inconsistency"], out["indirectness"], out["imprecision"], out["dependence"])
    if family == "effect" and worst > -2 and out["verdict"] == "Supported" and abs(out["est"]) >= THRESH["large"] \
            and min(abs(out["lo"]), abs(out["hi"])) >= THRESH["important"]:
        out["upgrade"] = 1
        reasons.append(f"upgrade: large effect (SMD {out['est']:.2f}, CI {out['lo']:.2f} to {out['hi']:.2f})")
    bg = BACKGROUND.get(c.claim_id)
    if bg and worst > -2 and out["verdict"].startswith("Supported"):
        out["upgrade"] = min(2, out["upgrade"] + 1)
        refs = "; ".join(r["citation"].split(".")[0] + " " + re.search(r"(\d{4})\.?$", r["citation"]).group(1)
                         for r in bg["references"])
        reasons.append(f"upgrade: triangulation with {bg['evidence_type']} ({bg['id']}: {refs})")
    out["not_assessed"] = NOT_ASSESSED
    total = out["rob"] + out["inconsistency"] + out["indirectness"] + out["imprecision"] + out["dependence"] + out["upgrade"]
    out["certainty"] = LEVELS[int(np.clip(level + total, 0, 3))]
    if out["verdict"] == "Untested":
        out["certainty"] = "Very low"
    out["reasons"] = "; ".join(reasons) if reasons else "no serious concerns in the domains assessed"
    # ---- scope of validity
    held = inf[inf.zo > 0] if expected != 0 else inf[inf.outcome.eq("reproduced")]
    out["scope_tested"] = ", ".join(sorted({x for v in inf.regions_tested.dropna() for x in str(v).split(";") if x}))
    out["scope_held"] = ", ".join(sorted({x for v in held.regions_tested.dropna() for x in str(v).split(";") if x}))
    flags_rows = inf.assign(claim_id=c.claim_id)[["claim_id", "dataset", "role", "template_used", "outcome", "z", "zo",
                                                  "n_case", "n_ref", "n_tests", "cohort", "flag", "weight"] +
                                                 (["d", "se_d", "conv"] if family == "effect" else [])]
    return out, [flags_rows]


BACKGROUND = {}  # claim id -> reference-knowledge record that directly supports it (claims/background.yaml)


def main():
    import argparse
    import yaml
    global DB, OUT
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=str(DB), help="claim store to rate (default from config/project.toml)")
    ap.add_argument("--out", default=str(OUT), help="folder for claim_rating.tsv and companions")
    a = ap.parse_args()
    DB, OUT = Path(a.db), Path(a.out)
    bgf = ROOT / "claims" / "background.yaml"
    for rec in (yaml.safe_load(bgf.read_text()) or []) if bgf.exists() else []:
        for cid in rec.get("supports_claims", []):
            BACKGROUND[cid] = rec
    claims, tests, ds, cov, old = load()
    cohort_of = ds.set_index("dataset").cohort.to_dict()
    tests = pick_variant(tests)
    by = dict(tuple(tests.groupby("claim_id")))
    rows, per_test = [], []
    for c in claims.itertuples():
        r, ft = rate_claim(c, by.get(c.claim_id, tests.iloc[:0]), cohort_of, cov)
        rows.append(r)
        per_test += ft
    rt = pd.DataFrame(rows)
    rt["old_tier"] = rt.claim_id.map(old)
    rt["rating"] = rt.verdict + ", " + rt.certainty.str.lower() + " certainty"
    rt["rating_version"] = RATING_VERSION
    OUT.mkdir(parents=True, exist_ok=True)
    rt.round(4).to_csv(OUT / "claim_rating.tsv", sep="\t", index=False)
    (pd.concat(per_test).round(4) if per_test else pd.DataFrame(columns=["claim_id", "dataset"])).to_csv(
        OUT / "rating_tests.tsv", sep="\t", index=False)
    xt = pd.crosstab(rt.old_tier, [rt.verdict, rt.certainty])
    xt.to_csv(OUT / "rating_vs_tiers.tsv", sep="\t")
    print(rt.verdict.value_counts().to_string(), "\n")
    print(pd.crosstab(rt.verdict, rt.certainty).reindex(columns=LEVELS[::-1], fill_value=0).to_string(), "\n")
    print(pd.crosstab(rt.old_tier, rt.verdict).to_string(), "\n")
    for d in ("rob", "inconsistency", "indirectness", "imprecision", "dependence", "upgrade"):
        print(d, rt[d].value_counts().sort_index().to_dict())


if __name__ == "__main__":
    main()

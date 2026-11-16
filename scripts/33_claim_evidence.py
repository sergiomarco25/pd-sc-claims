"""Evidence ranking of all published claims, with the context needed to read them (regions, donors, datasets).

Usage: python -I scripts/33_claim_evidence.py
Inputs: claims/claims.tsv, reports/R4_cross_claims.tsv (robustness per claim), results/tests/{internal,cross}/<ds>.tsv
(one row per claim x dataset: outcome, z, case and reference units, regions and classes tested),
"donors_*_tested" count the Lewy-body cases and controls of all datasets that could test the claim;
results/biology/claim_themes.tsv and donor_scores.tsv.gz (donors and regions per dataset).
Evidence tiers: established (robust), supported, unresolved (not replicated, inconsistent), contradicted, untested
(fewer than two other datasets could test it). Within a tier, claims are ranked by the number of datasets reproducing
them, then by the combined p-value. Writes results/meta/claim_evidence.tsv, claims/claim_evidence.tsv (git copy) and manuscript/table1_claims.md
(the 25 best-supported disease-effect claims as a markdown table).
"""
from pathlib import Path
import textwrap
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

RES = cfg.RESULTS
TIER = {"robust": "1 established", "supported": "2 supported", "not_replicated": "3 unresolved",
        "inconsistent": "3 unresolved", "contradicted": "4 contradicted", "insufficient": "5 untested"}
DISEASE_TYPES = {"composition", "de_gene", "de_pathway", "state_enrichment", "pathology_association",
                 "clinical_association", "genetic_enrichment"}


def tests(mode):
    t = pd.concat([pd.read_csv(f, sep="\t") for f in sorted((RES / "tests" / mode).glob("*.tsv"))], ignore_index=True)
    gw = RES / "tests" / f"{mode}_gwas"
    if gw.exists():
        g = pd.concat([pd.read_csv(f, sep="\t") for f in sorted(gw.glob("*.tsv"))], ignore_index=True)
        t = pd.concat([t[~t.set_index(["claim_id", "dataset"]).index.isin(g.set_index(["claim_id", "dataset"]).index)], g])
    return t


PAPER = {"smajic2022": "Smajić 2022", "kamath2022": "Kamath 2022", "lee2023": "Lee 2023", "martirosyan2024": "Martirosyan 2024",
         "wang2024": "Wang 2024", "adams2024": "Adams 2024", "ma2025": "Ma 2025", "volpato2025": "Volpato 2025",
         "menon2026": "Menon 2026", "barbareyes2026": "Barba-Reyes 2026", "barbareyes2025": "Barba-Reyes 2025",
         "xu2023": "Xu 2023", "prashant2024": "Prashant 2024", "feleke2021": "Feleke 2021", "zhu2024": "Zhu 2024",
         "chia2025": "Chia 2025", "jin2024": "Jin 2024", "glock2026": "Glock 2026", "lee2026": "Lee 2026",
         "clarence2025": "Clarence 2025", "dehestani2024": "Dehestani 2024"}
SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")


def fmt_p(x):
    if x == 0:
        return "<10⁻³⁰⁰"
    if x < 1e-3:
        m, ex = f"{x:.0e}".split("e")
        return f"{m}×10⁻{str(int(ex[1:])).translate(SUP)}"
    return f"{x:.3f}"


def table1(out, n=25):
    d = out[out.kind.eq("disease effect") & out.tier.str[0].isin(["1", "2"])].head(n)
    rows = ["| # | Claim | Paper | Claimed in | Reproduced / tested | Same region | Cases / controls tested | Combined p |",
            "|---|---|---|---|---|---|---|---|"]
    for i, r in enumerate(d.itertuples(), 1):
        sr = f"{int(r.same_region_rep)}/{int(r.same_region_k)}" if r.same_region_k and r.same_region_k > 0 else "–"
        rows.append(f"| {i} | {textwrap.shorten(r.statement, 120, placeholder='…')} | {PAPER.get(r.paper, r.paper)} | "
                    f"{r.region}, {r.cell_type} | {int(r.n_reproduced)}/{int(r.k)} | {sr} | "
                    f"{int(r.donors_cases_tested)}/{int(r.donors_controls_tested)} | {fmt_p(r.p)} |")
    (ROOT / "manuscript" / "table1_claims.md").write_text("\n".join(rows) + "\n")


def main():
    c = pd.read_csv(ROOT / "claims" / "claims.tsv", sep="\t").rename(columns={"id": "claim_id"})
    r = pd.read_csv(ROOT / "reports" / "R4_cross_claims.tsv", sep="\t")
    # claims that no dataset can test (spatial, clinical, chromatin-only ...) are absent from R4: add them as untested
    miss = c[~c.claim_id.isin(r.claim_id)]
    year = r.drop_duplicates("paper").set_index("paper").year
    r = pd.concat([r, miss[["claim_id", "paper", "statement", "centrality", "claim_type", "direction", "comparison",
                            "unit_of_analysis", "method", "quote", "location"]].assign(
        year=miss.paper.map(year), internal="not_testable", robustness="insufficient", k=0, n_reproduced=0,
        n_direction_only=0, n_not_reproduced=0, n_contradicted=0)], ignore_index=True)
    th = pd.read_csv(RES / "biology" / "claim_themes.tsv", sep="\t").groupby("claim_id").theme.agg("; ".join)
    s = pd.read_csv(RES / "biology" / "donor_scores.tsv.gz", sep="\t", usecols=["dataset", "region", "donor", "disease"])
    ds_info = s.drop_duplicates(["dataset", "donor"]).groupby("dataset").agg(
        donors=("donor", "nunique"), cases=("disease", lambda d: d.isin(["PD", "PDD", "DLB", "LBD"]).sum()),
        controls=("disease", lambda d: d.eq("control").sum()))
    ds_info["regions"] = s.groupby("dataset").region.agg(lambda x: ", ".join(sorted(set(x))))

    cross = tests("cross")
    cross = cross[~cross.outcome.isin(["not_testable"]) & cross.outcome.notna()]
    rep = cross[cross.outcome.eq("reproduced")]
    agg = cross.groupby("claim_id").agg(
        case_units=("n_case", "sum"), ref_units=("n_case", lambda x: cross.loc[x.index, "n_ref"].sum()),
        regions_tested=("regions_tested", lambda x: ", ".join(sorted({q for v in x.dropna() for q in str(v).split(";")}))),
        datasets_tested=("dataset", lambda x: ", ".join(sorted(x))))
    dn = cross.drop_duplicates(["claim_id", "dataset"])
    agg["donors_cases_tested"] = dn.groupby("claim_id").dataset.agg(lambda x: int(ds_info.cases.reindex(x).sum()))
    agg["donors_controls_tested"] = dn.groupby("claim_id").dataset.agg(lambda x: int(ds_info.controls.reindex(x).sum()))
    agg["datasets_reproducing"] = rep.groupby("claim_id").apply(
        lambda d: ", ".join(f"{a} ({b})" for a, b in zip(d.dataset, d.regions_tested.fillna("")))).reindex(agg.index)
    agg["datasets_contradicting"] = cross[cross.outcome.eq("contradicted")].groupby("claim_id").dataset.agg(", ".join).reindex(agg.index)

    internal = tests("internal").drop_duplicates("claim_id").set_index("claim_id")[["outcome", "z", "n_case", "n_ref"]]
    internal.columns = ["own_data_outcome", "own_data_z", "own_data_case_units", "own_data_ref_units"]

    out = (r.merge(c[["claim_id", "region", "cell_type", "cell_state", "genes", "reported_stat"]], on="claim_id", how="left")
           .merge(agg, left_on="claim_id", right_index=True, how="left")
           .merge(internal, left_on="claim_id", right_index=True, how="left"))
    for s in ["same_region", "other_region", "later_datasets"]:  # R4 stores the reproducing share: make it a count
        out[f"{s}_rep"] = (out[f"{s}_rep"] * out[f"{s}_k"]).round()
    out["themes"] = out.claim_id.map(th)
    dis = out.comparison.fillna("").str.contains(r"control|PD|LBD|DLB|Braak|stage|aged|disease", case=False)
    dis &= ~out.comparison.fillna("").str.contains(r"vs other|within", case=False)  # state-vs-rest contrasts in cases
    out["kind"] = np.where(out.claim_type.isin(DISEASE_TYPES) & dis, "disease effect", "descriptive")
    src = out.paper.map(lambda p: p)  # dataset id = paper id for papers with open data
    out["source_donors"] = src.map(ds_info.donors)
    out["source_cases"] = src.map(ds_info.cases)
    out["source_controls"] = src.map(ds_info.controls)
    out["source_regions"] = src.map(ds_info.regions)
    out["tier"] = out.robustness.map(TIER)
    out["share_reproducing"] = out.n_reproduced / out.k.replace(0, np.nan)
    out = out.sort_values(["tier", "n_reproduced", "p"], ascending=[True, False, True]).reset_index(drop=True)
    out.insert(0, "rank", np.arange(1, len(out) + 1))
    cols = ["rank", "tier", "kind", "claim_id", "paper", "year", "statement", "themes", "claim_type", "direction",
            "comparison", "region", "cell_type", "centrality", "unit_of_analysis", "reported_stat",
            "source_regions", "source_donors", "source_cases", "source_controls",
            "own_data_outcome", "own_data_z", "own_data_case_units", "own_data_ref_units",
            "k", "n_reproduced", "n_direction_only", "n_not_reproduced", "n_contradicted", "share_reproducing", "z", "p",
            "I2", "same_region_k", "same_region_rep", "other_region_k", "other_region_rep", "later_datasets_k",
            "later_datasets_rep", "donors_cases_tested", "donors_controls_tested", "case_units", "ref_units", "regions_tested", "datasets_tested",
            "datasets_reproducing", "datasets_contradicting", "robustness", "quote", "location"]
    out = out[cols]
    out.to_csv(RES / "meta" / "claim_evidence.tsv", sep="\t", index=False)
    out.to_csv(ROOT / "claims" / "claim_evidence.tsv", sep="\t", index=False)
    table1(out)
    print(out.groupby(["tier", "kind"]).size().to_string())
    pd.set_option("display.width", 250)
    pd.set_option("display.max_colwidth", 90)
    print(out[out.kind.eq("disease effect")].head(40)[["rank", "claim_id", "statement", "region", "k", "n_reproduced",
                                                       "donors_cases_tested", "donors_controls_tested"]].to_string(index=False))


if __name__ == "__main__":
    main()

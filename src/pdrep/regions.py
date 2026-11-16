"""Claim region vocabulary -> harmonised dataset regions (shared by claim contexts and cross-dataset tests)."""
REGION_MAP = {"SN": {"SN"}, "midbrain": {"midbrain", "SN"}, "striatum": {"putamen", "caudate", "striatum"},
              "putamen": {"putamen"}, "caudate": {"caudate"}, "DLPFC": {"DLPFC", "PFC"}, "PFC": {"PFC", "DLPFC"},
              "cingulate": {"cingulate"}, "temporal": {"temporal", "MTG"}, "MTG": {"MTG", "temporal"},
              "cerebellum": {"cerebellum"}}


def target_regions(claim_regions, have):
    """Regions of a dataset to test a claim in, and whether they match the claim's region ('same' / 'other')."""
    claim_regions = [r for r in claim_regions if r]
    if not claim_regions or "multiple" in claim_regions:
        return sorted(have), "same"
    want = set().union(*[REGION_MAP.get(r, set()) for r in claim_regions])
    hit = sorted(have & want)
    return (hit, "same") if hit else (sorted(have), "other")

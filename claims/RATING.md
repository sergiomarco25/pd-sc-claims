# Evidence rating rules, version 2.0.0

How every claim in the store gets a **verdict** (what the evidence says) and a **certainty** (how much to trust that
verdict), GRADE-style. This document is the citable statement of the rules (COMMUNITY_PLAN step 2.4); the
implementation is `scripts/52_rating.py`, whose `RATING_VERSION` must equal the version above and is written to
every row of `claims/store/claim_rating.tsv`. Design and motivation: `claims/ROADMAP.md` section 3.

Change policy: a change that can alter any verdict or certainty raises the minor version (2.1.0); a new threshold
value raises the patch version only if no rating changes; a change of the verdict or certainty categories raises the
major version. Each release of the store records the rating version it used.

## 1. Inputs

For each claim, the informative re-tests: rows of `claim_dataset` with outcome `reproduced`, `direction_only`,
`not_reproduced` or `contradicted` and a finite z. Per dataset one variant is used (preference: `both`, then
`revised`, then `preprint`). Tests that could not be run (`not_testable`) count for nothing. A claim without any
informative test is **Untested, very low certainty**: no dataset could test it, which is not the same as false.

Each test's z is **oriented**: positive means in the claimed direction (taken from the test engine's outcome);
direction-free claims keep the sign of the data.

## 2. Evidence families

| family | templates | how datasets are combined |
|---|---|---|
| effect | composition (CLR), DE genes, gene-set shift, state-gene shift, trajectory genes, communication genes | oriented z converted to a standardised mean difference, pooled by DerSimonian–Laird random effects |
| pattern | state existence, concordance between two diseases, GWAS (MAGMA) enrichment and specificity, familial-gene overlap | weighted Stouffer combination of oriented z, plus the share of datasets detecting the pattern |

A claim is in the pattern family only if all its tests use pattern templates.

**Effect size.** For a two-group test, d = z·√(1/n₁ + 1/n₂)/√m, where m is the number of region × cell-class tests
the engine combined; d is clipped to ±3 and SE = √(1/n₁ + 1/n₂ + d²/(2(n₁+n₂))). For a correlation with n units
(n ≥ 5), r = z/√n/√m (clipped to ±0.95), d = 2r/√(1−r²), SE = √(4/((n−1)(1−r²))). Findings that report a native
standardised effect and SE use those directly.

**Pattern weights.** w = √(n_case + n_ref) per test; combined z = Σ w·z / √Σ w².

## 3. Verdict

With the pooled SMD and its 95% CI [lo, hi] (effect family), equivalence margin 0.2:

| verdict | effect family | pattern family |
|---|---|---|
| Supported | CI excludes 0 on the claimed side (lo > 0) | detected (outcome reproduced) in ≥ 2/3 of datasets, or combined z > 1.96 with detection ≥ 1/2 |
| Partly supported | — | not Supported, but ≥ 2/3 of datasets point in the claimed direction |
| No effect | CI inside (−0.2, 0.2) | — |
| Contradicted | CI excludes 0 on the opposite side (hi < 0) | combined z < −1.96, or > 1/2 of datasets contradict |
| Inconclusive | anything else | anything else |
| Untested | no informative test | no informative test |

Claims without a direction: Supported if the CI excludes 0 on either side. **Null claims** (the claim asserts no
effect): "Supported (no effect)" if the CI lies inside (−0.2, 0.2) (pattern family: detected in ≥ 2/3 of datasets),
Contradicted if the CI excludes 0, else Inconclusive.

## 4. Certainty

**Start** from the number of independent cohorts with informative tests: ≥ 3 → High, 2 → Moderate, 1 → Low.
Then add the domain adjustments below (each 0, −1 or −2), add upgrades, and clip to the four levels
High / Moderate / Low / Very low. Untested claims are always Very low.

| domain | −1 (serious) | −2 (very serious) |
|---|---|---|
| Risk of bias | > 50% of weight from sorted nuclei in a composition test (datasets in `[rating] enriched_sort` of the project config) | — |
| Inconsistency (effect) | I² ≥ 50%, or < 70% of datasets in the pooled direction, or the 95% prediction interval crosses 0 for a Supported / Contradicted verdict | I² > 75% and < 60% of datasets in the pooled direction |
| Inconsistency (pattern) | detected in < 50% of datasets while the verdict is Supported or Partly supported | — |
| Indirectness | > 50% of weight from tests with one mismatch: other region than claimed, a partial proxy for the claimed cell state (signature coverage partial / text only / none), or sorted nuclei | > 50% of weight with two mismatches, or non-human evidence in the claim |
| Imprecision | effect: CI spans both trivial (< 0.2) and important (> 0.5) effects, or < 200 donor units; pattern: < 200 donor units | a single test, or < 30 donor units |
| Dependence | > 50% of weight from the source study's own cohort | all weight from the source cohort |

Agreement for inconsistency is judged against the **pooled** direction, so a contradicted claim whose datasets agree
on the opposite direction counts as consistent evidence.

**Upgrades** (only if no domain is −2): +1 for a large effect (Supported, |SMD| ≥ 0.8 and both CI limits ≥ 0.5 in
magnitude); +1 for triangulation with reference knowledge from another evidence type that directly supports the
claim (`claims/background.yaml`, `supports_claims`); at most +2 in total.

## 5. Scope and caveats

Every rating carries the regions tested and the regions where the evidence held (scope of validity), the reasons for
each adjustment in plain words, and caveats such as possible donor overlap between datasets drawn from the same
brain banks (reported, not merged).

Not yet assessed (listed on every rating, as worded in `[rating] not_assessed` of the project config; here: a gradient
with disease stage, triangulation beyond the reference records, and adjustment for age, sex, post-mortem interval and
RNA quality). The template families (section 2) are `[rating] effect_templates` and `pattern_templates`.

## 6. Thresholds

All thresholds are initial values for the sensitivity analysis, kept together in `THRESH` in `scripts/52_rating.py`:

| name | value | used for |
|---|---|---|
| equiv, trivial | 0.2 | equivalence margin (No effect), trivial effect |
| important | 0.5 | important effect (imprecision, upgrade) |
| large | 0.8 | large effect (upgrade) |
| i2_serious, i2_very | 0.5, 0.75 | inconsistency |
| conc_serious, conc_very | 0.7, 0.6 | share of datasets in the pooled direction |
| donors_ois, donors_min | 200, 30 | optimal information size, minimum donor units |
| weight_major | 0.5 | share of weight that triggers a domain |
| detect_supported | 2/3 | pattern family detection |

"""Claim re-tests (claimstore engine). One entry point, `test_claim(ctx, dataset, sigs)`, used for internal
replication (dataset = the claim's own) and for cross-dataset replication; `run(mode, dataset)` tests every claim of
claims/contexts.tsv in one dataset.

Inputs: a row of claims/contexts.tsv, a harmonised dataset, and that claim's rows of claims/signatures_all.tsv.gz.
Evidence comes from (a) the edgeR DE cache (results/de_cache/<dataset>/de/<region>__<class>__<contrast>.tsv.gz) for
gene-level claims, (b) cell metadata for composition, (c) pseudobulk signature scores for states and pathology.
Every result row carries an `outcome` from `classify_outcome`:
  reproduced      significant in the claimed direction (FDR < 0.05 for single genes / set p < 0.05)
  direction_only  claimed direction, not significant
  contradicted    significant in the opposite direction
  not_reproduced  opposite direction, not significant (or no expected direction and not significant)
  not_testable    data, groups or genes missing (reason in `note`)
Units of analysis are donor x region throughout. Multiple regions / classes are tested separately and combined
with a weighted Stouffer z (weights = sqrt(n units)).

Facts about the field come from the project's config: the stage ([stage] in config/contrasts.toml: its contrast code
for a continuous association and its label), and [tests] in config/project.toml:
  gene_aliases   gene symbols renamed since publication (MARCH*/SEPT* -> MARCHF*/SEPTIN* is built in)
  concordance    contrast codes that compare two contrasts' fold changes: code = {contrasts = [a, b], template = ".."}
  genetic_test   "module:function" testing genetic_enrichment claims (the field's genetics; none: not testable)
  unit_columns   cell metadata the composition test reads (default: the harmonised schema's donor columns)

Command line (analysis environment):
  python -I src/claimstore/cli.py tests internal DATASET   claims whose source dataset is DATASET
  python -I src/claimstore/cli.py tests cross DATASET      all other papers' claims, in DATASET
  add "gwas" to run only genetic_enrichment claims (written to <mode>_gwas/; "nogwas": all others)
Writes results/tests/<mode>/<dataset>.tsv.
"""
import importlib
import sys
from functools import lru_cache

import numpy as np
import pandas as pd
from scipy import stats

from claimstore import config as cfg
from claimstore.contrasts import STAGE, stage_value, unit_groups

ROOT = cfg.REPO
H = cfg.HARMONISED
RES = cfg.RESULTS
DE = RES / "de_cache"
UP = {"up", "enriched", "associated_pos", "exists"}
DOWN = {"down", "depleted", "associated_neg"}
TESTS = cfg.CONFIG.get("tests", {})
ALIASES = {str(k).upper(): str(v) for k, v in TESTS.get("gene_aliases", {}).items()}
CONCORDANCE = {code: (tuple(v["contrasts"]), v["template"]) for code, v in TESTS.get("concordance", {}).items()}
# [tests] signatures_per_class = true: a signature whose rows name cell classes is applied per class (each class gets
# the rows naming it, and rows naming none or "All"). Off by default: one gene->sign map for all of a claim's classes.
PER_CLASS = bool(TESTS.get("signatures_per_class", False))
STAGE_CONTRAST = str(STAGE.get("contrast", "stage"))
STAGE_LABEL = str(STAGE.get("label", "stage"))
N_PERM = 2000
RNG = np.random.default_rng(0)


def _schema_unit_columns():
    try:
        import tomllib
        with open(cfg.REPO / "config" / "harmonise.toml", "rb") as fh:
            spec = tomllib.load(fh)
    except OSError:
        spec = {}
    obs = spec.get("obs", [])
    donor = [c for c in obs if c not in ("dataset", "sample", "author_cell_type")]
    return list(dict.fromkeys(donor + ["enrichment"] + list(spec.get("donor_extras", []))))


UNIT_COLUMNS = list(TESTS.get("unit_columns", []) or _schema_unit_columns())


def sign_of(direction):
    return 1 if direction in UP else -1 if direction in DOWN else 0


def fname(cls):
    return cls.replace(" ", "_").replace("(", "").replace(")", "")


def classify_outcome(z, p, expected, alpha=0.05):
    if not np.isfinite(z) or not np.isfinite(p):
        return "not_testable"
    if expected == 0:
        return "reproduced" if p < alpha else "not_reproduced"
    s = np.sign(z) * expected
    if p < alpha:
        return "reproduced" if s > 0 else "contradicted"
    return "direction_only" if s > 0 else "not_reproduced"


def stouffer(zs, ns):
    zs, w = np.asarray(zs, float), np.sqrt(np.asarray(ns, float))
    ok = np.isfinite(zs)
    if not ok.any():
        return np.nan, np.nan
    z = (w[ok] * zs[ok]).sum() / np.sqrt((w[ok] ** 2).sum())
    return z, 2 * stats.norm.sf(abs(z))


def map_genes(genes, universe):
    up = {g.upper(): g for g in universe}
    out = []
    for g in genes:
        g0 = str(g).strip()
        for cand in (g0, ALIASES.get(g0.upper(), g0), g0.replace("MARCH", "MARCHF") if g0.startswith("MARCH") else g0,
                     g0.replace("SEPT", "SEPTIN") if g0.startswith("SEPT") and not g0.startswith("SEPTIN") else g0):
            if cand.upper() in up:
                out.append(up[cand.upper()])
                break
    return list(dict.fromkeys(out))


@lru_cache(maxsize=256)
def de_table(dataset, region, cls, contrast):
    f = DE / dataset / "de" / f"{region}__{fname(cls)}__{contrast}.tsv.gz"
    if not f.exists():
        return None
    t = pd.read_csv(f, sep="\t").set_index("gene")
    t["z"] = np.sign(t.logFC) * stats.norm.isf(np.clip(t.p, 1e-300, 1) / 2)
    return t


@lru_cache(maxsize=64)
def de_log(dataset):
    f = DE / dataset / "de_log.tsv"
    return pd.read_csv(f, sep="\t") if f.exists() else pd.DataFrame()


def gene_set_test(t, genes_dir):
    """genes_dir: dict gene -> expected sign (+1/-1/0). Returns dict of statistics on one DE table."""
    g = [x for x in genes_dir if x in t.index]
    if not g:
        return None
    e = np.array([genes_dir[x] for x in g], float)
    e_eff = np.where(e == 0, 1.0, e)  # no expected direction: test for up
    z = t.loc[g, "z"].values * e_eff
    obs = z.mean()
    if len(g) == 1:
        zz = t.loc[g[0], "z"] * e_eff[0]
        return {"n_genes": 1, "z": zz, "p": float(t.loc[g[0], "padj"]), "concordance": float(zz > 0),
                "frac_sig_concordant": float((t.loc[g[0], "padj"] < 0.05) and zz > 0),
                "effect": float(t.loc[g[0], "logFC"] * e_eff[0]), "test": "single gene (FDR)"}
    allz = t.z.values
    null = np.array([allz[RNG.choice(len(allz), len(g), replace=False)].mean() for _ in range(N_PERM)])
    # random sets ignore the signs, so compare |set mean| with a sign-symmetric null
    zz = (obs - null.mean()) / null.std(ddof=1)
    p = 2 * stats.norm.sf(abs(zz))
    sig = t.loc[g, "padj"].values < 0.05
    return {"n_genes": len(g), "z": zz, "p": p, "concordance": float((z > 0).mean()),
            "frac_sig_concordant": float((sig & (z > 0)).mean()), "effect": float(np.median(t.loc[g, "logFC"].values * e_eff)),
            "test": "competitive set test (random sets)"}


@lru_cache(maxsize=32)
def cell_obs(dataset):
    import anndata as ad
    obs = ad.read_h5ad(H / f"{dataset}.h5ad", backed="r").obs
    cols = [c for c in UNIT_COLUMNS if c in obs]
    o = obs[cols].copy()
    if "enrichment" in o:
        o = o[~o.enrichment.astype(str).str.lower().eq("sorted")]
    return o


def unit_meta(dataset):
    o = cell_obs(dataset)
    o = o.assign(unit=o.donor.astype(str) + "|" + o.region.astype(str))
    meta = o.drop(columns=["cell_class"]).groupby("unit", observed=True).first()
    return meta.join(unit_groups(meta))


def composition_test(dataset, regions, classes, contrast, expected):
    o = cell_obs(dataset)
    o = o.assign(unit=o.donor.astype(str) + "|" + o.region.astype(str))
    meta = unit_meta(dataset)
    counts = pd.crosstab(o.unit, o.cell_class.astype(str))
    prop = (counts + 0.5).div((counts + 0.5).sum(1), axis=0)
    clr = np.log(prop).sub(np.log(prop).mean(1), axis=0)
    rows = []
    for reg in regions:
        for cls in classes:
            if cls not in clr:
                continue
            u = meta.index[meta.region.astype(str).eq(reg)]
            if contrast == STAGE_CONTRAST:
                b = stage_value(meta.loc[u]).dropna()
                if len(b) < 6 or b.nunique() < 3:
                    continue
                rho, p = stats.spearmanr(clr.loc[b.index, cls], b)
                rows.append({"region": reg, "cell_class": cls, "z": np.sign(rho) * stats.norm.isf(p / 2), "p": p,
                             "effect": rho, "n_case": len(b), "n_ref": 0})
                continue
            g = meta.loc[u, contrast] if contrast in meta else pd.Series("", index=u)
            a, r = clr.loc[g.index[g.eq("case")], cls], clr.loc[g.index[g.eq("ref")], cls]
            if len(a) < 3 or len(r) < 3:
                continue
            p = stats.mannwhitneyu(a, r, alternative="two-sided").pvalue
            d = a.mean() - r.mean()
            rows.append({"region": reg, "cell_class": cls, "z": np.sign(d) * stats.norm.isf(p / 2), "p": p,
                         "effect": d, "n_case": len(a), "n_ref": len(r),
                         "prop_case": float(prop.loc[a.index, cls].mean()), "prop_ref": float(prop.loc[r.index, cls].mean())})
    return pd.DataFrame(rows)


MIN_SCORE_GENES = 5


@lru_cache(maxsize=64)
def unit_table(dataset):
    f = DE / dataset / "input" / "units.tsv"
    return pd.read_csv(f, sep="\t", dtype=str).set_index("unit") if f.exists() else None


@lru_cache(maxsize=32)
def unit_logcpm_z(dataset, region, cls):
    """Pseudobulk log2(CPM+1) of one region x class, z-scored per gene across units (genes x units)."""
    f = DE / dataset / "input" / f"{region}__{fname(cls)}.csv.gz"
    if not f.exists():
        return None
    x = pd.read_csv(f, index_col=0)
    lc = np.log2(x.div(x.sum(0), axis=1) * 1e6 + 1)
    return lc.sub(lc.mean(1), axis=0).div(lc.std(1).replace(0, np.nan), axis=0).dropna()


def donor_score_test(dataset, region, cls, contrast, genes_dir):
    """Self-contained set test with donors as the unit: per-unit score = mean over the set's genes of
    (expected sign x z-scored log CPM); case vs reference by Mann-Whitney (or Spearman with the stage).
    Positive z = as claimed. Respects between-gene correlation, unlike gene-permutation set tests."""
    zt, ut = unit_logcpm_z(dataset, region, cls), unit_table(dataset)
    if zt is None or ut is None:
        return None
    hit = {}
    for g, sgn in genes_dir.items():
        m = map_genes([g], zt.index)
        if m:
            hit[m[0]] = sgn if sgn != 0 else 1
    if len(hit) < MIN_SCORE_GENES:
        return None
    score = zt.loc[list(hit)].mul(pd.Series(hit), axis=0).mean(0)
    if contrast == STAGE_CONTRAST:
        b = stage_value(ut.loc[score.index]).dropna()
        if len(b) < 6 or b.nunique() < 3:
            return None
        rho, p = stats.spearmanr(score.loc[b.index], b)
        return {"z": np.sign(rho) * stats.norm.isf(max(p, 1e-300) / 2), "p": p, "effect": rho, "n_case": len(b),
                "n_ref": 0, "n_genes": len(hit), "test": f"donor score ~ {STAGE_LABEL} (Spearman)"}
    if contrast not in ut:
        return None
    g = ut.loc[score.index, contrast].fillna("")
    a, r = score[g.eq("case").values], score[g.eq("ref").values]
    if len(a) < 3 or len(r) < 3:
        return None
    p = stats.mannwhitneyu(a, r, alternative="two-sided").pvalue
    d = a.mean() - r.mean()
    return {"z": np.sign(d) * stats.norm.isf(max(p, 1e-300) / 2), "p": p, "effect": d, "n_case": len(a),
            "n_ref": len(r), "n_genes": len(hit), "test": "donor signature score (Mann-Whitney)"}


def gene_level_test(dataset, regions, classes, contrast, genes_dir):
    """Per region x class: donor-level signature-score test for sets of >= MIN_SCORE_GENES genes (primary);
    edgeR FDR for single genes; gene-permutation set test only for 2-4 genes. Concordance and share of genes
    significant in the claimed direction are always reported from the edgeR table."""
    rows = []
    for reg in regions:
        for cls in classes:
            t = de_table(dataset, reg, cls, contrast)
            if t is None:
                continue
            mapped = {}
            for g, sgn in genes_dir.items():  # map gene by gene so each keeps its own expected sign
                m = map_genes([g], t.index)
                if m:
                    mapped[m[0]] = sgn
            r = gene_set_test(t, mapped)
            if r is None:
                continue
            ds = donor_score_test(dataset, reg, cls, contrast, mapped) if len(mapped) >= MIN_SCORE_GENES else None
            if ds is not None:
                r.update({k: ds[k] for k in ["z", "p", "effect", "test"]})
                r["n_case"], r["n_ref"] = ds["n_case"], ds["n_ref"]
            else:
                lg = de_log(dataset)
                n = lg[(lg.region == reg) & (lg.cell_class == cls) & (lg.contrast == contrast)]
                r["n_case"] = int(n.n_case.iloc[0]) if len(n) else np.nan
                r["n_ref"] = int(n.n_ref.iloc[0]) if len(n) and pd.notna(n.n_ref.iloc[0]) else 0
            r.update({"region": reg, "cell_class": cls})
            rows.append(r)
    return pd.DataFrame(rows)


def gene_level_by_class(dataset, regions, classes, contrast, sig, claim_sign, per_class=None):
    """gene_level_test with one gene -> sign map for every class, or (per_class, default PER_CLASS) one map per class
    from the signature rows that name it; rows naming no class (or "All") go to every class."""
    per_class = PER_CLASS if per_class is None else per_class
    named = sig.get("cell_type", pd.Series("", index=sig.index)).fillna("").astype(str).str.split(";").map(
        lambda xs: {x.strip() for x in xs if x.strip() and x.strip() != "All"})
    if not per_class or "ALL" in classes or not named.map(bool).any():
        return gene_level_test(dataset, regions, classes, contrast, expected_signs(sig, claim_sign))
    parts = []
    for cls in classes:
        rows = sig[named.map(lambda c: not c or cls in c)]
        if len(rows):
            parts.append(gene_level_test(dataset, regions, [cls], contrast, expected_signs(rows, claim_sign)))
    return pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()


def expected_signs(sig, claim_sign):
    out = {}
    for g, k in zip(sig.gene, sig.kind):
        out[g] = 1 if k == "de_up" else -1 if k == "de_down" else claim_sign
    return out


def summarise(df, expected, label, extra=None):
    base = {"template_used": label, **(extra or {})}
    if df is None or df.empty:
        return {**base, "outcome": "not_testable", "note": "no region/class/contrast with enough units or genes"}
    n = df.n_case.fillna(0) + df.n_ref.fillna(0)
    if len(df) == 1:
        z, p = float(df.z.iloc[0]), float(df.p.iloc[0])
    else:
        z, p = stouffer(df.z, n)
    out = {**base, "z": z, "p": p, "effect": float(np.nanmedian(df.effect)), "n_tests": len(df),
           "set_test": ";".join(sorted(df.test.dropna().unique())) if "test" in df else "",
           "regions_tested": ";".join(sorted(df.region.astype(str).unique())),
           "classes_tested": ";".join(sorted(df.cell_class.astype(str).unique())),
           "n_case": int(df.n_case.max()) if df.n_case.notna().any() else 0,
           "n_ref": int(df.n_ref.max()) if df.n_ref.notna().any() else 0}
    for c in ["n_genes", "concordance", "frac_sig_concordant", "prop_case", "prop_ref"]:
        if c in df:
            out[c] = float(np.nanmean(df[c]))
    out["outcome"] = classify_outcome(z, p, expected)
    return out


MAX_STATE_CELLS = 30_000
N_RANDOM_SETS = 200


@lru_cache(maxsize=8)
def class_cells(dataset, regions, classes):
    """Log-normalised expression of up to MAX_STATE_CELLS unsorted nuclei of the given classes/regions (tuples)."""
    import anndata as ad
    import scanpy as sc
    a = ad.read_h5ad(H / f"{dataset}.h5ad", backed="r")
    obs = a.obs
    m = obs.region.astype(str).isin(regions).values
    if "ALL" not in classes:
        m &= obs.cell_class.astype(str).isin(classes).values
    if "enrichment" in obs:
        m &= ~obs.enrichment.astype(str).str.lower().eq("sorted").values
    idx = np.where(m)[0]
    if len(idx) > MAX_STATE_CELLS:
        idx = np.sort(RNG.choice(idx, MAX_STATE_CELLS, replace=False))
    if len(idx) < 200:
        return None
    b = a[idx].to_memory()
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    return b


@lru_cache(maxsize=32)
def class_specificity(dataset):
    f = RES / "gwas" / "sets" / f"{dataset}.spec.tsv.gz"
    return pd.read_csv(f, sep="\t", index_col=0) if f.exists() else None


def state_existence_test(dataset, regions, classes, genes, max_genes=50):
    """Do the state's markers co-vary within the cell class beyond random genes matched for mean expression AND
    specificity to the parent class, after removing each gene's dependence on nuclear library size? And do the
    high-scoring nuclei come from many donors? (Depth and lineage genes alone cannot make a state pass.)"""
    b = class_cells(dataset, tuple(regions), tuple(classes))
    if b is None:
        return None
    X = b.X.tocsc() if hasattr(b.X, "tocsc") else b.X
    names = np.asarray(b.var_names)
    pos = {g: i for i, g in enumerate(names)}
    mean = np.asarray(X.mean(0)).ravel()
    det = np.asarray((X > 0).mean(0)).ravel()
    sel = [g for g in map_genes(genes, names) if det[pos[g]] >= 0.01][:max_genes]
    gi = [pos[g] for g in sel]
    if len(gi) < 3:
        return None
    depth = np.log1p(np.asarray(b.obs["total_counts"] if "total_counts" in b.obs else X.sum(1)).ravel().astype(float))
    depth = (depth - depth.mean()) / (depth.std() + 1e-9)

    def resid(ix):
        M = X[:, ix].toarray() if hasattr(X, "toarray") else np.asarray(X[:, ix])
        M = M - M.mean(0)
        return M - np.outer(depth, depth @ M / len(depth))

    def coherence(ix):
        C = np.corrcoef(resid(ix), rowvar=False)
        return np.nanmean(C[np.triu_indices(len(ix), 1)])

    spec = class_specificity(dataset)
    cls = [c for c in classes if c != "ALL"]
    sp = np.zeros(len(names))
    if spec is not None and cls and all(c in spec.columns for c in cls):
        sp = spec.reindex(names)[cls].max(1).fillna(0).values
    eb = pd.qcut(pd.Series(mean).rank(method="first"), 10, labels=False).values
    sb = pd.qcut(pd.Series(sp).rank(method="first"), 5, labels=False).values if sp.any() else np.zeros(len(names), int)
    key = eb * 10 + sb
    ok = det >= 0.01
    pools = {k: np.where((key == k) & ok)[0] for k in np.unique(key)}
    obs_c = coherence(gi)
    null = np.array([coherence([int(RNG.choice(pools[key[i]])) if len(pools[key[i]]) else i for i in gi])
                     for _ in range(N_RANDOM_SETS)])
    z = (obs_c - null.mean()) / null.std(ddof=1)
    p = (np.sum(null >= obs_c) + 1) / (len(null) + 1)
    R = resid(gi)
    score = (R / (R.std(0) + 1e-9)).mean(1)
    # discreteness: share of nuclei above the 99th percentile of scores from random matched sets (expected 1%)
    null_scores = []
    for _ in range(20):
        Rn = resid([int(RNG.choice(pools[key[i]])) if len(pools[key[i]]) else i for i in gi])
        null_scores.append((Rn / (Rn.std(0) + 1e-9)).mean(1))
    thr = np.quantile(np.concatenate(null_scores), 0.99)
    tail_ratio = float((score > thr).mean() / 0.01)
    top = score >= np.quantile(score, 0.9)
    share = pd.Series(top).groupby(b.obs.donor.astype(str).values).mean()
    spread = float((share >= 0.02).mean())
    return {"region": ";".join(regions), "cell_class": ";".join(classes), "z": z, "p": p, "effect": obs_c - null.mean(),
            "n_genes": len(gi), "donor_spread": spread, "tail_ratio": tail_ratio, "n_case": int(len(share)), "n_ref": 0,
            "test": "state co-expression and high-score tail vs expression+specificity-matched random genes, depth-adjusted"}


def concordance_test(dataset, regions, classes, c1, c2):
    """Spearman correlation of the fold changes of two contrasts, per region x class (a [tests] concordance code)."""
    rows = []
    for reg in regions:
        for cls in classes:
            t1, t2 = de_table(dataset, reg, cls, c1), de_table(dataset, reg, cls, c2)
            if t1 is None or t2 is None:
                continue
            g = t1.index.intersection(t2.index)
            if len(g) < 100:
                continue
            rho, p = stats.spearmanr(t1.loc[g, "logFC"], t2.loc[g, "logFC"])
            rows.append({"region": reg, "cell_class": cls, "z": np.sign(rho) * stats.norm.isf(max(p, 1e-300) / 2),
                         "p": p, "effect": rho, "n_case": len(g), "n_ref": 0})
    return pd.DataFrame(rows)


def genetic_hook():
    """The field's test of genetic_enrichment claims ([tests] genetic_test = "module:function"), or None."""
    spec = str(TESTS.get("genetic_test", "")).strip()
    if not spec:
        return None
    mod, fn = spec.split(":")
    sys.path.insert(0, str(cfg.REPO / "src"))
    return getattr(importlib.import_module(mod), fn)


def genetic_test(ctx, dataset, statement, modality, cell_state):
    fn = genetic_hook()
    if fn is None:
        return [{"template_used": "genetic_enrichment", "outcome": "not_testable",
                 "note": "no genetic test for this field ([tests] genetic_test in config/project.toml)"}]
    return fn(ctx, dataset, statement, modality, cell_state)


def test_claim(ctx, dataset, sigs, claim=None):
    """Returns a list of result dicts (one per signature variant). `claim`: the claims.yaml row (for GWAS claims).

    "No difference" claims (contexts.null_claim): reproduced if the test is not significant (p >= 0.05), contradicted
    if it is (equivalence-style rule; lenient where data are underpowered, as reported)."""
    res = _test_claim(ctx, dataset, sigs, claim)
    if str(getattr(ctx, "null_claim", False)).lower() == "true":
        for r in res:
            if r.get("outcome") not in (None, "not_testable") and pd.notna(r.get("p")):
                r["outcome"] = "reproduced" if r["p"] >= 0.05 else "contradicted"
                r["note"] = (str(r.get("note") or "") + " null claim: equivalence-style rule").strip()
    return res


def _test_claim(ctx, dataset, sigs, claim=None):
    regions = [r for r in str(ctx.regions).split(";") if r]
    classes = [c for c in str(ctx.cell_classes).split(";") if c]
    exp = sign_of(ctx.direction)
    con = ctx.contrast
    if ctx.test_template == "genetic_enrichment" and claim is not None:
        return genetic_test(ctx, dataset, str(claim.statement), claim.modality, ctx.cell_state)
    if con in CONCORDANCE:
        (c1, c2), label = CONCORDANCE[con]
        df = concordance_test(dataset, regions, [c for c in classes if c != "ALL"] or ["ALL"], c1, c2)
        return [summarise(df, exp if exp else 1, label)]
    if con == "within" and ctx.test_template in ("state_signature", "geneset_shift", "de_pseudobulk", "composition_donor"):
        if sigs is None or sigs.empty:
            return [{"template_used": "state_existence", "outcome": "not_testable", "note": "no signature genes"}]
        out = []
        for variant, sv in sigs.groupby(sigs.get("variant", pd.Series("both", index=sigs.index)).fillna("both")):
            pick = sv[sv.kind.isin(["state_markers", "named_genes", "de_up"])]
            r = state_existence_test(dataset, regions, [c for c in classes] or ["ALL"], list(pick.gene.unique()))
            if r is None:
                out.append({"template_used": "state_existence", "variant": variant, "outcome": "not_testable",
                            "note": "too few nuclei or marker genes detected"})
                continue
            res = summarise(pd.DataFrame([r]), 1, "state_existence", {"variant": variant})
            res["donor_spread"], res["tail_ratio"] = r["donor_spread"], r["tail_ratio"]
            if res["outcome"] == "reproduced" and r["tail_ratio"] < 2:
                res["outcome"], res["note"] = "direction_only", "coherent programme, no distinct subpopulation"
            elif res["outcome"] == "reproduced" and r["donor_spread"] < 0.5:
                res["outcome"], res["note"] = "direction_only", "distinct but concentrated in few donors"
            out.append(res)
        return out
    if con in ("within", "sex"):
        return [{"template_used": ctx.test_template, "outcome": "not_testable",
                 "note": f"contrast '{con}' with template '{ctx.test_template}' handled in a later pass"}]
    if ctx.test_template == "genetic_enrichment":
        if claim is None:
            return [{"template_used": "gwas_magma", "outcome": "not_testable", "note": "claim text not passed"}]
        return genetic_test(ctx, dataset, str(claim.statement), claim.modality, ctx.cell_state)
    cons = con.split("+")
    if ctx.test_template == "composition_donor" and not str(ctx.cell_state).strip():
        df = pd.concat([composition_test(dataset, regions, classes, c, exp) for c in cons], ignore_index=True)
        return [summarise(df, exp, "composition_clr")]
    if ctx.test_template in ("de_pseudobulk", "geneset_shift", "state_signature", "composition_donor",
                             "pathology_correlation", "communication_liana", "trajectory"):
        # communication / trajectory claims with a group contrast: expression-based proxy on their genes
        # (ligands and receptors, or trajectory genes, in the claim's cell types)
        if sigs is None or sigs.empty:
            return [{"template_used": ctx.test_template, "outcome": "not_testable", "note": "no signature genes"}]
        out = []
        for variant, sv in sigs.groupby(sigs.get("variant", pd.Series("both", index=sigs.index)).fillna("both")):
            cls = classes
            if "ALL" in cls or not cls:
                cls = ["ALL"]
            signs = expected_signs(sv, exp)
            df = pd.concat([gene_level_by_class(dataset, regions, cls, c, sv, exp) for c in cons], ignore_index=True)
            # gene z-scores are already multiplied by each gene's expected sign: positive = as claimed
            oriented = 1 if any(v != 0 for v in signs.values()) else 0
            label = {"de_pseudobulk": "de_genes", "geneset_shift": "geneset_shift",
                     "communication_liana": "communication_genes", "trajectory": "trajectory_genes"}.get(
                ctx.test_template, "state_genes_shift")
            out.append(summarise(df, oriented, label, {"variant": variant,
                                                 "signatures": ";".join(sorted(sv.signature_id.unique()))[:400]}))
        return out
    return [{"template_used": ctx.test_template, "outcome": "not_testable",
             "note": f"template '{ctx.test_template}' handled in a later pass"}]


SIG_COLUMNS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "gene", "gene_reported", "variant"]


def signatures():
    """claims/signatures_all.tsv.gz, one row per claim x gene (empty when the project has no signatures yet)."""
    f = ROOT / "claims" / "signatures_all.tsv.gz"
    sig = (pd.read_csv(f, sep="\t", dtype=str, keep_default_na=False) if f.exists()
           else pd.DataFrame(columns=SIG_COLUMNS, dtype=str))
    sig = sig.assign(claim_id=sig.claim_ids.str.split(";")).explode("claim_id")
    sig["claim_id"] = sig.claim_id.str.strip()
    sig["gene"] = sig.gene.where(sig.gene != "", sig.gene_reported)
    return sig


def run(mode, ds, part="all"):
    """Test claims in one dataset: mode "internal" (claims whose source dataset is ds) or "cross" (all other claims,
    regions retargeted to ds); part "gwas" / "nogwas" / "all" (genetic_enrichment claims only / all but them / all).
    Writes results/tests/<mode>[_gwas]/<ds>.tsv and returns the table."""
    import yaml
    from claimstore.contexts import target_regions
    ctx = pd.read_csv(ROOT / "claims" / "contexts.tsv", sep="\t").fillna("")
    sig = signatures()
    ctx = ctx[ctx.test_template.eq("genetic_enrichment")] if part == "gwas" else \
        ctx[~ctx.test_template.eq("genetic_enrichment")] if part == "nogwas" else ctx
    claims_all = pd.DataFrame(yaml.safe_load((ROOT / "claims" / "claims.yaml").read_text())).set_index("id")
    todo = ctx[ctx.dataset.eq(ds)] if mode == "internal" else ctx[~ctx.dataset.eq(ds)]
    if mode == "cross":  # retarget regions to this dataset; skip claims that are not testable anywhere
        have = set(pd.read_csv(H / f"{ds}_donors.tsv", sep="\t").region.astype(str))
        # a claim whose own dataset is not harmonised here ("no dataset") is still tested in every other dataset
        todo = todo[todo.status.isin(["ok", "no dataset"])].assign(status="ok")
        reg = todo.claim_id.map(lambda i: target_regions(
            claims_all.region.get(i) if isinstance(claims_all.region.get(i), list) else [claims_all.region.get(i)], have))
        todo["regions"], todo["region_match"] = reg.str[0].str.join(";"), reg.str[1]
    rows = []
    for c in todo.itertuples():
        if c.status != "ok":
            rows.append({"claim_id": c.claim_id, "dataset": ds, "outcome": "not_testable", "note": c.status})
            continue
        extra = {"region_match": c.region_match} if mode == "cross" else {}
        for r in test_claim(c, ds, sig[sig.claim_id.eq(c.claim_id)], claims_all.loc[c.claim_id]):
            rows.append({"claim_id": c.claim_id, "dataset": ds, "contrast": c.contrast, **extra, **r})
    out = RES / "tests" / (mode + ("_gwas" if part == "gwas" else ""))
    out.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows, columns=None if rows else ["claim_id", "dataset", "outcome"])
    df.to_csv(out / f"{ds}.tsv", sep="\t", index=False)
    print(ds, mode, len(df), "rows")
    if len(df):
        print(df.outcome.value_counts().to_string())
    return df


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if len(argv) not in (2, 3) or argv[0] not in ("internal", "cross") or (len(argv) == 3 and
                                                                            argv[2] not in ("gwas", "nogwas", "all")):
        print("usage: tests internal|cross DATASET [gwas|nogwas]")
        return 2
    run(argv[0], argv[1], argv[2] if len(argv) == 3 else "all")
    return 0


if __name__ == "__main__":
    sys.exit(main())

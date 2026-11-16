"""Run the analysis a question needs, from a declarative spec, on every harmonised dataset (question loop, roadmap step 3).

Usage: python -I scripts/54_question_analysis.py <question_id>
Reads the spec from claims/questions.yaml (`proposed.spec`), writes results/questions/<question_id>.tsv in the format
pdrep.store.resolve_evidence() reads (one row per dataset x region x cell class: est, se, n_case, n_ref or n).
The same scoring as scripts/27_programmes.py: per donor x region pseudobulk (results/de_cache/<ds>/input), log2 CPM,
z-scored across units, programme score = mean z of its genes; groups compared by Hedges' g, covariates by Spearman
rho (Fisher z).

Spec (YAML), one of `genes`, `gene_set` or `signature`, and one of `compare` or `correlate`:
  genes: [CIRBP, RBM3]                         explicit genes
  gene_set: {library: h.all | c2.cp.reactome | c5.go.bp, name: HALLMARK_HYPOXIA}
  signature: barbareyes2026__astro_UPR_module1 a signature id of claims/signatures_all.tsv.gz
  cell_classes: [Oligodendrocytes]            harmonised classes (ALL = all nuclei)
  regions: [SN, midbrain]                      optional; default all regions
  compare: {contrast: PDany_vs_control}        an existing contrast column of units.tsv, or
  compare: {column: sex, case: [M], ref: [F], subset: {disease: [PD, PDD, DLB, LBD]}}
  correlate: {column: age, subset: {disease: [control]}}   numeric donor column (age, pmi, braak_lb, unified_lb_stage)
  min_genes: 5
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

BASE = cfg.WORK
DE, GS, OUT = BASE / "results" / "de_cache", BASE / "raw" / "genesets", BASE / "results" / "questions"


def genes_of(spec):
    if spec.get("genes"):
        return [g.upper() for g in spec["genes"]]
    if spec.get("gene_set"):
        lib, name = spec["gene_set"]["library"], spec["gene_set"]["name"]
        for line in open(GS / f"{lib}.v2024.1.Hs.symbols.gmt"):
            f = line.rstrip("\n").split("\t")
            if f[0] == name:
                return [g.upper() for g in f[2:]]
        raise KeyError(f"{name} not in {lib}")
    if spec.get("signature"):
        sig = pd.read_csv(ROOT / "claims" / "signatures_all.tsv.gz", sep="\t", dtype=str, usecols=["signature_id", "gene"])
        return sig[sig.signature_id == spec["signature"]].gene.str.upper().unique().tolist()
    raise ValueError("spec needs genes, gene_set or signature")


def hedges(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    n1, n2 = len(a), len(b)
    if n1 < 3 or n2 < 3:
        return np.nan, np.nan
    sp = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1)) / (n1 + n2 - 2))
    if sp == 0:
        return np.nan, np.nan
    g = (1 - 3 / (4 * (n1 + n2) - 9)) * (a.mean() - b.mean()) / sp
    return g, np.sqrt((n1 + n2) / (n1 * n2) + g ** 2 / (2 * (n1 + n2)))


def fisher(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 6 or np.unique(y[ok]).size < 3:
        return np.nan, np.nan, int(ok.sum())
    rho = stats.spearmanr(x[ok], y[ok])[0]
    return np.arctanh(np.clip(rho, -0.999, 0.999)), 1 / np.sqrt(ok.sum() - 3), int(ok.sum())


def subset_mask(u, subset):
    m = pd.Series(True, index=u.index)
    for col, vals in (subset or {}).items():
        m &= u[col].astype(str).isin([str(v) for v in (vals if isinstance(vals, list) else [vals])])
    return m.values


def run(spec):
    genes = genes_of(spec)
    classes = [c.replace(" ", "_").replace("(", "").replace(")", "") for c in spec["cell_classes"]]
    regions = set(spec.get("regions") or [])
    min_genes = spec.get("min_genes", 2 if len(genes) <= 3 else 5)
    rows = []
    for ds_dir in sorted(DE.iterdir()):
        uf = ds_dir / "input" / "units.tsv"
        if not uf.exists():
            continue
        units = pd.read_csv(uf, sep="\t", dtype=str).set_index("unit")
        for f in sorted((ds_dir / "input").glob("*__*.csv.gz")):
            region, cls_f = f.name[: -len(".csv.gz")].split("__", 1)
            if cls_f not in classes or (regions and region not in regions):
                continue
            x = pd.read_csv(f, index_col=0)
            x.index = x.index.str.upper()
            lc = np.log2(x.div(x.sum(0), axis=1) * 1e6 + 1)
            z = lc.sub(lc.mean(1), axis=0).div(lc.std(1).replace(0, np.nan), axis=0)
            g = [x_ for x_ in genes if x_ in z.index]
            if len(g) < min_genes:
                continue
            score = z.loc[g].mean(0).values
            u = units.loc[z.columns]
            sub = subset_mask(u, (spec.get("compare") or spec.get("correlate") or {}).get("subset"))
            base = dict(dataset=ds_dir.name, region=region, cell_class=cls_f.replace("_", " "), n_genes=len(g))
            if "compare" in spec:
                c = spec["compare"]
                if "contrast" in c:
                    grp = u[c["contrast"]].fillna("").values if c["contrast"] in u else np.array([""] * len(u))
                    case, ref = (grp == "case") & sub, (grp == "ref") & sub
                else:
                    col = u[c["column"]].astype(str).values
                    case = np.isin(col, [str(v) for v in c["case"]]) & sub
                    ref = np.isin(col, [str(v) for v in c["ref"]]) & sub
                est, se = hedges(score[case], score[ref])
                if np.isfinite(est):
                    rows.append({**base, "test": c.get("contrast") or f"{c['column']}: {c['case']} vs {c['ref']}",
                                 "est": est, "se": se, "n_case": int(case.sum()), "n_ref": int(ref.sum())})
            else:
                c = spec["correlate"]
                cov = pd.to_numeric(u[c["column"]], errors="coerce").values.astype(float)
                cov[~sub] = np.nan
                est, se, n = fisher(score, cov)
                if np.isfinite(est):
                    rows.append({**base, "test": f"correlation with {c['column']}", "est": est, "se": se, "n": n})
    return pd.DataFrame(rows)


def main():
    qid = sys.argv[1]
    qs = yaml.safe_load((ROOT / "claims" / "questions.yaml").read_text()) or []
    q = next(x for x in qs if x["id"] == qid)
    spec = q["proposed"]["spec"]
    res = run(spec)
    OUT.mkdir(parents=True, exist_ok=True)
    res.to_csv(OUT / f"{qid}.tsv", sep="\t", index=False)
    print(f"{qid}: {len(res)} dataset x region x class units, {res.dataset.nunique() if len(res) else 0} datasets"
          f" -> {OUT / (qid + '.tsv')}")
    if len(res):
        print(res.groupby("dataset")[["est"]].mean().round(2).T.to_string())


if __name__ == "__main__":
    main()

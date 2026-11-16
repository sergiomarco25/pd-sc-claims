"""Structure of the disease programmes across patients, cell types and regions (donor-level, all datasets).

Usage: python -I scripts/30_donor_analyses.py
Input: results/biology/donor_scores.tsv.gz and donor_composition.tsv.gz (scripts/29_donor_scores.py). Cases are
donors with Lewy-body disease (PD, PDD, DLB, LBD); iLBD and other diagnoses are left out unless stated.
  A. heterogeneity   share of cases above the control 95th percentile (5% expected) and log SD ratio cases/controls
  B. coordination    Spearman correlation of a programme between cell types across donor units of one dataset x region
                     (all donors, cases only, controls only); share of variance on the first shared axis vs shuffled
  C. stress vs reactive  correlation of stress programmes with published reactive / inflammatory glial states (cases, controls)
  D. regions         in multi-region datasets: correlation of a donor's score between regions, and case-control effect
                     per region
  E. links           within cases: programme score vs Lewy stage, microglia share, excitatory and total neuron share
  F. stress index    per donor: mean over glial classes (astrocytes, microglia, oligodendrocytes, OPCs) and regions of
                     the heat-response and chaperone-folding scores, z-scored against the controls of each dataset x
                     region x class; compared between male and female cases (Hedges' g, pooled) and summarised by diagnosis
Effects are pooled across datasets by DerSimonian-Laird (correlations as Fisher z). Writes results/biology/donor_*.tsv.
"""
from pathlib import Path
import itertools
import warnings
import sys
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.meta import fixed, dersimonian_laird  # noqa: E402

from pdrep import config as cfg  # noqa: E402

B = cfg.RESULTS / "biology"
CASES = {"PD", "PDD", "DLB", "LBD"}
CORE = ["heat response", "chaperone folding", "unfolded protein response", "hypoxia", "cold-inducible (CIRBP, RBM3)",
        "UPR astrocytes (Barba-Reyes)", "UPR microglia (Barba-Reyes)", "reactive astrocytes (Barba-Reyes)",
        "reactive microglia (Barba-Reyes)", "GPNMB-high microglia (Smajić)", "DAM microglia (Barba-Reyes)", "AR2 astrocytes (Barba-Reyes)",
        "inflammatory response", "oxidative phosphorylation", "HSP90AA1-high microglia (Smajić)",
        "C3/CD44 astrocytes, Astrocytes1 (Martirosyan)", "UPR astrocytes, Astrocytes3 (Martirosyan)",
        "GFAP/MT astrocytes, Astrocytes4 (Martirosyan)", "stress oligodendrocytes, Oligos5 (Martirosyan)",
        "activated astrocytes, Ast-1 (Xu)", "activated microglia, Micr-1 (Xu)", "activated myeloid cells (Ma)",
        "disease-associated oligodendrocytes (Adams)"]
GLIA = ["Astrocytes", "Microglia", "Oligodendrocytes", "OPCs", "Endothelial", "Excitatory neurons", "Inhibitory neurons"]
RNG = np.random.default_rng(0)
warnings.simplefilter("ignore", FutureWarning)


def corr_fz(x, y, min_n=8):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < min_n:
        return np.nan, np.nan, int(ok.sum())
    r = stats.spearmanr(x[ok], y[ok])[0]
    return np.arctanh(np.clip(r, -0.999, 0.999)), 1 / np.sqrt(ok.sum() - 3), int(ok.sum())


def pool(df, keys):
    """DL pooling of est/se over datasets (dataset x region rows combined within dataset first)."""
    out = []
    for k, d in df.groupby(keys):
        per = d.groupby("dataset").apply(lambda x: pd.Series(fixed(x.est, x.se), index=["est", "se"]))
        r = dersimonian_laird(per.est, per.se)
        out.append({**dict(zip(keys, k if isinstance(k, tuple) else (k,))), **r, "n_units": int(d.n.sum()),
                    "datasets": ";".join(per.dropna().index)})
    return pd.DataFrame(out)


def heterogeneity(s):
    """Tails: share of cases beyond the control 95th (above) / 5th (below) percentile, against 5% by chance and against
    the share expected if all cases were shifted by the mean difference with the control spread (a pure shift)."""
    rows = []
    for (ds, reg, cls, pg), d in s.groupby(["dataset", "region", "cell_class", "programme"]):
        cd = d[d.disease.isin(CASES)]
        a, r = cd.score.values, d[d.disease.eq("control")].score.values
        if len(a) < 4 or len(r) < 4:
            continue
        hi, lo = np.quantile(r, 0.95), np.quantile(r, 0.05)
        mu, sd = r.mean(), r.std(ddof=1)
        sd_ratio = np.log(a.std(ddof=1) / sd) if sd > 0 else np.nan
        se = np.sqrt(1 / (2 * (len(a) - 1)) + 1 / (2 * (len(r) - 1)))
        g = (a.mean() - mu) / sd if sd > 0 else np.nan
        row = {"dataset": ds, "region": reg, "cell_class": cls, "programme": pg, "n_cases": len(a), "n_ctrl": len(r),
               "cases_above_ctrl95": int((a > hi).sum()), "cases_below_ctrl5": int((a < lo).sum()),
               "exp_above_shift": len(a) * stats.norm.sf((hi - mu) / sd - g) if sd > 0 else np.nan,
               "exp_below_shift": len(a) * stats.norm.cdf((lo - mu) / sd - g) if sd > 0 else np.nan,
               "est": sd_ratio, "se": se, "n": len(a) + len(r)}
        for dx in ["PD", "PDD", "DLB"]:
            x = cd.score.values[cd.disease.eq(dx).values]
            row[f"{dx}_above"], row[f"{dx}_n"] = int((x > hi).sum()), len(x)
        rows.append(row)
    h = pd.DataFrame(rows)
    p = pool(h, ["programme", "cell_class"]).rename(columns={"est": "log_sd_ratio"})
    cols = ["cases_above_ctrl95", "cases_below_ctrl5", "exp_above_shift", "exp_below_shift", "n_cases",
            "PD_above", "PD_n", "PDD_above", "PDD_n", "DLB_above", "DLB_n"]
    tot = h.groupby(["programme", "cell_class"])[cols].sum()
    tot["share_above"] = tot.cases_above_ctrl95 / tot.n_cases
    tot["share_below"] = tot.cases_below_ctrl5 / tot.n_cases
    tot["share_above_if_shift"] = tot.exp_above_shift / tot.n_cases
    tot["share_below_if_shift"] = tot.exp_below_shift / tot.n_cases
    tot["binom_p"] = [stats.binomtest(int(a), int(n), 0.05, alternative="greater").pvalue for a, n in zip(tot.cases_above_ctrl95, tot.n_cases)]
    # excess over a pure shift: Poisson-binomial approximated by Poisson with the expected count
    tot["excess_p"] = stats.poisson.sf(tot.cases_above_ctrl95 - 1, tot.exp_above_shift)
    return h, p.merge(tot.reset_index(), on=["programme", "cell_class"])


def coordination(s):
    rows, pc_rows = [], []
    for pg in ["heat response", "chaperone folding", "hypoxia", "cold-inducible (CIRBP, RBM3)", "oxidative phosphorylation"]:
        d = s[(s.programme == pg) & s.cell_class.isin(GLIA)]
        for (ds, reg), dd in d.groupby(["dataset", "region"]):
            w = dd.pivot_table(index="unit", columns="cell_class", values="score")
            dis = dd.drop_duplicates("unit").set_index("unit").disease.reindex(w.index)
            for grp, mask in [("all donors", np.ones(len(w), bool)), ("cases", dis.isin(CASES).values),
                              ("controls", dis.eq("control").values)]:
                for c1, c2 in itertools.combinations(sorted(w.columns), 2):
                    est, se, n = corr_fz(w[c1].values[mask], w[c2].values[mask])
                    if np.isfinite(est):
                        rows.append({"dataset": ds, "region": reg, "programme": pg, "group": grp, "pair": f"{c1} ~ {c2}",
                                     "est": est, "se": se, "n": n})
            m = w.dropna(axis=1, thresh=int(0.8 * len(w))).dropna()
            if m.shape[0] >= 8 and m.shape[1] >= 3:
                z = (m - m.mean()) / m.std()
                ev = np.linalg.svd(z.values, compute_uv=False) ** 2
                obs = ev[0] / ev.sum()
                null = []
                for _ in range(200):
                    zz = np.column_stack([RNG.permutation(z.values[:, j]) for j in range(z.shape[1])])
                    e = np.linalg.svd(zz, compute_uv=False) ** 2
                    null.append(e[0] / e.sum())
                pc_rows.append({"dataset": ds, "region": reg, "programme": pg, "n_units": m.shape[0], "n_classes": m.shape[1],
                                "pc1_share": obs, "null_share": float(np.mean(null)), "p": float((np.sum(np.array(null) >= obs) + 1) / 201)})
    c = pd.DataFrame(rows)
    return c, pool(c, ["programme", "group", "pair"]), pd.DataFrame(pc_rows)


STATES = [("heat response", "Astrocytes"), ("heat response", "Microglia"), ("heat response", "Oligodendrocytes"),
          ("UPR astrocytes (Barba-Reyes)", "Astrocytes"), ("UPR astrocytes, Astrocytes3 (Martirosyan)", "Astrocytes"),
          ("UPR microglia (Barba-Reyes)", "Microglia"), ("HSP90AA1-high microglia (Smajić)", "Microglia"),
          ("GPNMB-high microglia (Smajić)", "Microglia"), ("activated microglia, Micr-1 (Xu)", "Microglia"),
          ("activated myeloid cells (Ma)", "Microglia"), ("reactive microglia (Barba-Reyes)", "Microglia"),
          ("inflammatory response", "Microglia"), ("C3/CD44 astrocytes, Astrocytes1 (Martirosyan)", "Astrocytes"),
          ("activated astrocytes, Ast-1 (Xu)", "Astrocytes"), ("AR2 astrocytes (Barba-Reyes)", "Astrocytes"),
          ("reactive astrocytes (Barba-Reyes)", "Astrocytes"), ("inflammatory response", "Astrocytes")]


def upr_vs_reactive(s):
    """Correlation across donors between stress programmes and published reactive / inflammatory glial states of
    several studies (all pairs of STATES), in cases and controls separately."""
    rows = []
    for (p1, c1), (p2, c2) in itertools.combinations(STATES, 2):
        a = s[(s.programme == p1) & (s.cell_class == c1)][["dataset", "region", "unit", "disease", "score"]]
        b = s[(s.programme == p2) & (s.cell_class == c2)][["dataset", "region", "unit", "score"]]
        m = a.merge(b, on=["dataset", "region", "unit"], suffixes=("_1", "_2"))
        for grp, mm in [("cases", m[m.disease.isin(CASES)]), ("controls", m[m.disease.eq("control")])]:
            for (ds, reg), d in mm.groupby(["dataset", "region"]):
                est, se, n = corr_fz(d.score_1.values, d.score_2.values)
                if np.isfinite(est):
                    rows.append({"dataset": ds, "region": reg, "pair": f"{p1} [{c1}] ~ {p2} [{c2}]", "group": grp,
                                 "est": est, "se": se, "n": n})
    u = pd.DataFrame(rows)
    return u, pool(u, ["pair", "group"])


def regions(s):
    multi = s.groupby("dataset").region.nunique()
    multi = multi[multi >= 2].index
    rows, eff = [], []
    for pg in ["heat response", "chaperone folding", "UPR microglia (Barba-Reyes)", "cold-inducible (CIRBP, RBM3)", "hypoxia"]:
        d = s[(s.programme == pg) & s.dataset.isin(multi)]
        for (ds, cls), dd in d.groupby(["dataset", "cell_class"]):
            w = dd.pivot_table(index="donor", columns="region", values="score")
            dis = dd.drop_duplicates("donor").set_index("donor").disease.reindex(w.index)
            for r1, r2 in itertools.combinations(sorted(w.columns), 2):
                for grp, mask in [("all donors", np.ones(len(w), bool)), ("cases", dis.isin(CASES).values)]:
                    est, se, n = corr_fz(w[r1].values[mask], w[r2].values[mask])
                    if np.isfinite(est):
                        rows.append({"dataset": ds, "programme": pg, "cell_class": cls, "group": grp, "pair": f"{r1} ~ {r2}",
                                     "est": est, "se": se, "n": n})
        for (ds, reg, cls), dd in s[(s.programme == pg)].groupby(["dataset", "region", "cell_class"]):
            a, r = dd[dd.disease.isin(CASES)].score.values, dd[dd.disease.eq("control")].score.values
            if len(a) >= 3 and len(r) >= 3:
                sp = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(r) - 1) * r.var(ddof=1)) / (len(a) + len(r) - 2))
                g = (a.mean() - r.mean()) / sp if sp > 0 else np.nan
                eff.append({"dataset": ds, "region": reg, "cell_class": cls, "programme": pg, "hedges_g": g,
                            "se": np.sqrt((len(a) + len(r)) / (len(a) * len(r)) + g ** 2 / (2 * (len(a) + len(r)))),
                            "n_cases": len(a), "n_ctrl": len(r)})
    rc = pd.DataFrame(rows)
    return rc, pool(rc, ["programme", "cell_class", "group"]), pd.DataFrame(eff)


def links(s, comp):
    comp = comp.pivot_table(index=["dataset", "unit"], columns="cell_class", values="share")
    rows = []
    for pg, cls in [("heat response", "Oligodendrocytes"), ("heat response", "ALL"), ("chaperone folding", "Oligodendrocytes"),
                    ("UPR microglia (Barba-Reyes)", "Microglia"), ("UPR astrocytes (Barba-Reyes)", "Astrocytes"),
                    ("reactive microglia (Barba-Reyes)", "Microglia"), ("cold-inducible (CIRBP, RBM3)", "ALL"),
                    ("GPNMB-high microglia (Smajić)", "Microglia"), ("hypoxia", "Microglia"),
                    ("oxidative phosphorylation", "Excitatory neurons"), ("HSP90AA1-high microglia (Smajić)", "Microglia"),
                    ("activated myeloid cells (Ma)", "Microglia"), ("activated microglia, Micr-1 (Xu)", "Microglia"),
                    ("C3/CD44 astrocytes, Astrocytes1 (Martirosyan)", "Astrocytes")]:
        d = s[(s.programme == pg) & (s.cell_class == cls) & s.disease.isin(CASES)].copy()
        d = d.join(comp, on=["dataset", "unit"])
        for (ds, reg), dd in d.groupby(["dataset", "region"]):
            for var in ["stage", "age", "pmi", "Microglia", "Excitatory neurons", "neurons_total", "Astrocytes", "Oligodendrocytes"]:
                if var not in dd:
                    continue
                est, se, n = corr_fz(dd.score.values, pd.to_numeric(dd[var], errors="coerce").values, min_n=6)
                if np.isfinite(est):
                    rows.append({"dataset": ds, "region": reg, "programme": pg, "cell_class": cls, "variable": var,
                                 "est": est, "se": se, "n": n})
    lk = pd.DataFrame(rows)
    return lk, pool(lk, ["programme", "cell_class", "variable"])


def stress_index(s):
    x = s[s.programme.isin(["heat response", "chaperone folding"])
          & s.cell_class.isin(["Astrocytes", "Microglia", "Oligodendrocytes", "OPCs"])].copy()

    def zc(d):
        c = d[d.disease.eq("control")].score
        return (d.score - c.mean()) / c.std(ddof=1) if len(c) >= 4 else d.score * np.nan
    x["z"] = x.groupby(["dataset", "region", "cell_class", "programme"], group_keys=False).apply(zc)
    idx = x.groupby(["dataset", "donor", "disease", "sex"]).z.mean().dropna().reset_index()  # regions averaged per donor
    idx["high"] = idx.z > stats.norm.isf(0.05)
    cases = idx[idx.disease.isin(CASES)].copy()
    cases["sex"] = cases.sex.astype(str).str.upper().str[0]
    rows = []
    for ds, d in cases.groupby("dataset"):
        m, f = d[d.sex.eq("M")].z, d[d.sex.eq("F")].z
        if len(m) < 3 or len(f) < 3:
            continue
        sp = np.sqrt(((len(m) - 1) * m.var() + (len(f) - 1) * f.var()) / (len(m) + len(f) - 2))
        g = (m.mean() - f.mean()) / sp
        rows.append({"dataset": ds, "est": g, "se": np.sqrt(1 / len(m) + 1 / len(f) + g ** 2 / (2 * (len(m) + len(f)))),
                     "n": len(m) + len(f), "n_male": len(m), "n_female": len(f),
                     "high_male": d[d.sex.eq("M")].high.mean(), "high_female": d[d.sex.eq("F")].high.mean()})
    sx = pd.DataFrame(rows)
    pooled = pool(sx.assign(k="male vs female cases"), ["k"])
    summ = idx.assign(group=np.where(idx.disease.isin(CASES), idx.disease, idx.disease)).groupby("group").agg(
        donors=("z", "size"), mean_z=("z", "mean"), share_high=("high", "mean")).reset_index()
    for sxk in ["M", "F"]:
        d = cases[cases.sex.eq(sxk)]
        summ.loc[len(summ)] = {"group": f"cases {sxk}", "donors": len(d), "mean_z": d.z.mean(), "share_high": d.high.mean()}
    return idx, sx, pooled, summ


def main():
    s = pd.read_csv(B / "donor_scores.tsv.gz", sep="\t", low_memory=False)
    s = s[s.programme.isin(CORE)]
    comp = pd.read_csv(B / "donor_composition.tsv.gz", sep="\t")
    h, hp = heterogeneity(s)
    hp.to_csv(B / "donor_heterogeneity.tsv", sep="\t", index=False)
    h.to_csv(B / "donor_heterogeneity_per_dataset.tsv", sep="\t", index=False)
    c, cp, pc = coordination(s)
    cp.to_csv(B / "donor_coordination.tsv", sep="\t", index=False)
    pc.to_csv(B / "donor_coordination_pc1.tsv", sep="\t", index=False)
    u, up = upr_vs_reactive(s)
    up.to_csv(B / "donor_upr_vs_reactive.tsv", sep="\t", index=False)
    rc, rp, eff = regions(s)
    rp.to_csv(B / "donor_regions_corr.tsv", sep="\t", index=False)
    eff.to_csv(B / "donor_regions_effects.tsv", sep="\t", index=False)
    lk, lp = links(s, comp)
    lp.to_csv(B / "donor_links.tsv", sep="\t", index=False)
    idx, sx, sxp, summ = stress_index(s)
    idx.to_csv(B / "donor_stress_index.tsv", sep="\t", index=False)
    sx.to_csv(B / "donor_stress_sex.tsv", sep="\t", index=False)
    pd.concat([sxp, summ]).to_csv(B / "donor_stress_summary.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    show = lambda d, cols: print(d[cols].round(3).to_string(index=False), "\n")
    print("A. heterogeneity (glia, key programmes)")
    show(hp[hp.cell_class.isin(["Oligodendrocytes", "Astrocytes", "Microglia", "ALL"]) & hp.programme.isin(CORE[:9])].sort_values("binom_p"),
         ["programme", "cell_class", "k", "log_sd_ratio", "p", "share_above", "share_above_if_shift", "excess_p",
          "share_below", "share_below_if_shift", "PD_above", "PD_n", "PDD_above", "PDD_n", "DLB_above", "DLB_n"])
    print("B. coordination (heat response, pooled r = tanh(est))")
    x = cp[cp.programme.eq("heat response")].assign(r=lambda d: np.tanh(d.est))
    show(x.sort_values(["group", "r"], ascending=[True, False]), ["group", "pair", "k", "r", "p"])
    print("B2. first shared axis vs shuffled")
    show(pc.groupby("programme")[["pc1_share", "null_share"]].median().reset_index().assign(
        datasets_p05=pc.groupby("programme").p.apply(lambda p: int((p < 0.05).sum())).values,
        units=pc.groupby("programme").size().values), ["programme", "pc1_share", "null_share", "datasets_p05", "units"])
    print("C. UPR vs reactive")
    show(up.assign(r=np.tanh(up.est)), ["pair", "group", "k", "r", "p"])
    print("D. cross-region within donor")
    show(rp.assign(r=np.tanh(rp.est)), ["programme", "cell_class", "group", "k", "r", "p"])
    print("E. links within cases")
    show(lp.assign(r=np.tanh(lp.est)).sort_values("p").head(30), ["programme", "cell_class", "variable", "k", "r", "p"])
    print("F. donor stress index")
    print(sx.round(3).to_string(index=False))
    print(sxp.round(4).to_string(index=False))
    print(summ.round(3).to_string(index=False))


if __name__ == "__main__":
    main()

"""Where does the PD stress state sit on the AD trajectories of Green et al. 2024 (ROSMAP DLPFC, 437 donors)?

Usage: python -I scripts/50a_green_rosmap.py
Uses only the open supplementary tables of Green et al. 2024 (Nature 633:634; raw/green2024): per-donor snRNA-seq
state proportions (Supp. Table 3), clinicopathological traits (Supp. Table 1), Palantir trajectory probabilities and
pseudotime (Supp. Table 5) and the fitted state dynamics along each trajectory (Supp. Table 5).
Heat-shock states of the paper (HSPA1A/HSPH1/DNAJB1-high: Mic.11, Ast.9, End.5, Oli.8; Exc.1 adds immediate-early
genes) are compared with the PD stress programme; stress states without HSPs (Ast.10, Oli.7) and the lipid-associated
microglia (Mic.12, Mic.13) are references. Per state, sqrt proportion is regressed (OLS, standardised) on each trait with
age, sex and post-mortem interval as covariates; coordination is the Spearman correlation of heat-shock state
proportions across cell types, against all pairs of states from different cell types. A donor heat-shock index is the
mean z(sqrt proportion) of the four heat-shock states. Writes results/biology/green2024_rosmap.tsv and
green2024_coordination.tsv.
"""
from pathlib import Path
import itertools
import warnings
import numpy as np
import pandas as pd
from scipy import stats

warnings.simplefilter("ignore")

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

RAW = cfg.RAW / "green2024"
B = cfg.RESULTS / "biology"
HS = ["Mic.11", "Ast.9", "End.5", "Oli.8"]
SHOW = HS + ["Exc.1", "Ast.10", "Oli.7", "Mic.12", "Mic.13", "Ast.5", "Mic.10", "OPC.3", "Mic.6", "Mic.15", "Ast.1",
                   "Mic.14", "Ast.7"]
TRAITS = ["sqrt.amyloid_mf", "sqrt.tangles_mf", "cogng_demog_slope", "prAD", "ABA", "pseudotime", "pmi", "age_death"]


def load():
    x = lambda n, s: pd.read_excel(RAW / f"41586_2024_7871_MOESM{n}_ESM.xlsx", sheet_name=s)
    prop = x(5, "snRNA-seq proportions").set_index("individualID")
    tr = x(3, "Discovery cohort").set_index("individualID")
    pal = x(7, "Palantir Trajectories").set_index("individualID")[["pseudotime", "prAD", "ABA"]]
    dyn = x(7, "Dynamics predicted values")
    return prop, tr.join(pal, how="inner"), dyn


def ols(y, x, cov):
    d = pd.concat([y.rename("y"), x.rename("x"), cov], axis=1).dropna()
    d = (d - d.mean()) / d.std()
    X = np.c_[np.ones(len(d)), d.drop(columns="y").values]
    b, *_ = np.linalg.lstsq(X, d.y.values, rcond=None)
    r = d.y.values - X @ b
    dof = len(d) - X.shape[1]
    se = np.sqrt(np.diag(r @ r / dof * np.linalg.inv(X.T @ X)))
    return b[1], se[1], 2 * stats.t.sf(abs(b[1] / se[1]), dof), len(d)


def main():
    prop, tr, dyn = load()
    sq = np.sqrt(prop.reindex(tr.index))
    z = (sq - sq.mean()) / sq.std()
    sq["heat-shock index"] = z[HS].mean(1)
    rows = []
    for s in SHOW + ["heat-shock index"]:
        for t in TRAITS:
            cov = tr[[c for c in ["age_death", "msex", "pmi"] if c != t]]
            b, se, p, n = ols(sq[s], tr[t], cov)  # state ~ trait + covariates
            rows.append({"state": s, "trait": t, "beta": b, "se": se, "p": p, "n": n})
    # dynamics: change from start to end of each trajectory (fitted proportion, log2 ratio)
    for s in SHOW:
        for traj, d in dyn[dyn.feature == s].groupby("trajectory"):
            d = d.sort_values("pseudotime")
            a, e = d.fit.iloc[:3].mean(), d.fit.iloc[-3:].mean()
            rows.append({"state": s, "trait": f"dynamics {traj} (log2 end/start)", "beta": np.log2(max(e, 1e-4) / max(a, 1e-4))})
    t = pd.DataFrame(rows)
    t["fdr"] = np.nan
    ok = t.p.notna()
    t.loc[ok, "fdr"] = stats.false_discovery_control(t.loc[ok, "p"])
    t.to_csv(B / "green2024_rosmap.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(t.pivot_table(index="state", columns="trait", values="beta").reindex(SHOW + ["heat-shock index"]).round(2).to_string())
    print(t[t.p.notna()].pivot_table(index="state", columns="trait", values="fdr").reindex(SHOW + ["heat-shock index"]).applymap(lambda v: f"{v:.1e}").to_string())
    # coordination of heat-shock states across cell types
    r = sq[prop.columns].corr("spearman")
    ct = lambda s: s.split(".")[0]
    pairs = [(a, b) for a, b in itertools.combinations(prop.columns, 2) if ct(a) != ct(b)]
    allr = np.array([r.loc[a, b] for a, b in pairs])
    hs = [(a, b, r.loc[a, b]) for a, b in itertools.combinations(HS, 2)]
    c = pd.DataFrame(hs, columns=["a", "b", "rho"])
    c["percentile_among_cross_type_pairs"] = [(allr < v).mean() * 100 for v in c.rho]
    c.to_csv(B / "green2024_coordination.tsv", sep="\t", index=False)
    print(c.round(3).to_string(index=False))
    print("median |rho| all cross-type pairs", np.median(np.abs(allr)).round(3))
    # donors by terminal fate: heat-shock index in late prAD vs late ABA vs early
    late = tr.pseudotime > tr.pseudotime.quantile(0.66)
    grp = np.where(~late, "early", np.where(tr.prAD > tr.ABA, "late prAD", "late ABA"))
    print(sq.assign(g=grp).groupby("g")[["heat-shock index"] + HS].mean().round(3).to_string())
    for a, b in [("late prAD", "early"), ("late ABA", "early"), ("late prAD", "late ABA")]:
        print(a, "vs", b, stats.mannwhitneyu(sq["heat-shock index"][grp == a], sq["heat-shock index"][grp == b]).pvalue)


if __name__ == "__main__":
    main()

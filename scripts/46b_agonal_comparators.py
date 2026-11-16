"""Do non-neurodegenerative deaths with terminal hypoxia or fever produce the same stress state?

Usage: python -I scripts/46b_agonal_comparators.py
Two open comparators without Lewy pathology:
  - Fullard et al. 2021 (Genome Med; GSE164485): COVID-19 (5) vs controls (4); medulla oblongata, prefrontal cortex,
    choroid plexus; cells x genes counts (dense CSV) with author cell types.
  - Wei et al. 2023 (Nat Commun; GSE240457): opioid overdose vs drug-free controls; ventral midbrain (substantia
    nigra + VTA); pseudobulk per donor x cell type extracted by scripts/46a_wei2023_extract.R (used when present).
Per dataset x region x cell type, donor pseudobulks are scored as in scripts/27_programmes.py (mean z-scored log2 CPM
of the programme's genes across donors) and cases are compared with controls (Hedges' g, Mann-Whitney p).
For Wei et al., programme scores are also correlated with brain pH across donors (an agonal marker; lower pH after
prolonged agonal states). Writes results/biology/agonal_comparators.tsv.
"""
from pathlib import Path
import importlib.util
import sys
import warnings
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
warnings.simplefilter("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

RAW = cfg.RAW
B = cfg.RESULTS / "biology"
spec = importlib.util.spec_from_file_location("prog", ROOT / "scripts" / "27_programmes.py")
prog = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prog)
HEAT10 = ["HSPA1A", "HSPA1B", "HSPH1", "DNAJB1", "DNAJA1", "HSP90AA1", "BAG3", "CHORDC1", "STIP1", "AHSA1"]
SHOW = ["shared heat-shock genes (10)", "heat response", "chaperone folding", "unfolded protein response", "hypoxia",
        "cold-inducible (CIRBP, RBM3)", "disease-associated oligodendrocytes (Adams)", "UPR microglia (Barba-Reyes)",
        "UPR astrocytes (Barba-Reyes)", "inflammatory response"]


def fullard():
    meta = pd.read_csv(RAW / "fullard2021covid" / "GSE164485_meta_data_20210104.csv.gz", index_col=0)
    key = meta.subject_id.astype(str) + "|" + meta.tissue.astype(str) + "|" + meta.celltype.astype(str)
    sums = None
    for ch in pd.read_csv(RAW / "fullard2021covid" / "GSE164485_gene_expression_20210104.csv.gz", index_col=0,
                          chunksize=5000):
        part = ch.groupby(key.reindex(ch.index).values).sum()
        sums = part if sums is None else sums.add(part, fill_value=0)
        print("fullard chunk", flush=True)
    pb = sums.T
    pb.index = pb.index.str.upper()
    cols = pd.DataFrame([c.split("|") for c in pb.columns], columns=["donor", "region", "cell_class"], index=pb.columns)
    cov = meta.drop_duplicates("subject_id").set_index("subject_id").is_covid.astype(int)
    cols["case"] = cols.donor.map(cov).astype(bool)
    nuc = key.value_counts()
    cols["nuclei"] = nuc.reindex(pb.columns).values
    return pb, cols.assign(dataset="fullard2021 (COVID-19)")


def wei():
    f = RAW / "wei2023" / "wei2023_pseudobulk.csv.gz"
    if not f.exists():
        return None, None
    pb = pd.read_csv(f, index_col=0)
    pb.index = pb.index.str.upper()
    don = pd.read_csv(RAW / "wei2023" / "wei2023_donors.tsv", sep="\t").set_index("sampleID")
    cols = pd.DataFrame([c.split("|") for c in pb.columns], columns=["donor", "cell_class"], index=pb.columns)
    cols["case"] = cols.donor.map(don.category).eq("Opioid")
    cols["ph"] = cols.donor.map(don.ph).astype(float)
    cols["region"] = "ventral midbrain"
    cols["nuclei"] = np.nan
    return pb, cols.assign(dataset="wei2023 (opioid overdose)")


def score(pb, cols, min_donors=3):
    progs = {k: v[0] for k, v in prog.programmes().items()}
    progs["shared heat-shock genes (10)"] = HEAT10
    rows = []
    for (ds, reg, cls), c in cols.groupby(["dataset", "region", "cell_class"]):
        x = pb[c.index]
        x = x.loc[:, x.sum(0) > 2e4]
        c = c.loc[x.columns]
        if c.case.sum() < min_donors or (~c.case).sum() < min_donors:
            continue
        lc = np.log2(x.div(x.sum(0), axis=1) * 1e6 + 1)
        lc = lc[(x > 0).mean(1) >= 0.5]
        z = lc.sub(lc.mean(1), axis=0).div(lc.std(1).replace(0, np.nan), axis=0)
        for p in SHOW:
            g = [k for k in progs[p] if k in z.index]
            if len(g) < (2 if len(progs[p]) == 2 else 5):
                continue
            s = z.loc[g].mean(0)
            a, r = s[c.case.values], s[~c.case.values]
            sp = np.sqrt(((len(a) - 1) * a.var() + (len(r) - 1) * r.var()) / (len(a) + len(r) - 2))
            row = {"dataset": ds, "region": reg, "cell_class": cls, "programme": p, "n_case": len(a),
                   "n_ctrl": len(r), "hedges_g": (a.mean() - r.mean()) / sp if sp > 0 else np.nan,
                   "p": stats.mannwhitneyu(a, r).pvalue}
            if "ph" in c and c.ph.notna().sum() >= 10:  # brain pH, an agonal marker: correlation across all donors
                ok = c.ph.notna().values
                rho, pr = stats.spearmanr(s.values[ok], c.ph.values[ok])
                row.update({"rho_ph": rho, "p_ph": pr})
            rows.append(row)
    return pd.DataFrame(rows)


def main():
    out = []
    for f in (fullard, wei):
        pb, cols = f()
        if pb is not None:
            out.append(score(pb, cols))
    t = pd.concat(out, ignore_index=True)
    t.to_csv(B / "agonal_comparators.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    keep = t.cell_class.str.contains("Oli|Ast|Mic|olig|astro|micro|OPC|Endo", case=False)
    print(t[keep].pivot_table(index=["dataset", "region", "cell_class"], columns="programme", values="hedges_g")
          .reindex(columns=SHOW).round(2).to_string())


if __name__ == "__main__":
    main()

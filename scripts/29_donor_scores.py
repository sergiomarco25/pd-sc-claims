"""Donor-level programme scores: one row per dataset x donor x region x cell class x programme.

Usage: python -I scripts/29_donor_scores.py
Same scoring as scripts/27_programmes.py (mean z-scored log2 CPM of a programme's genes across the donor units of one
dataset x region x cell class pseudobulk), but the per-unit scores are kept, together with donor metadata, so that
patient heterogeneity, coordination between cell types, cross-region consistency and links to composition can be
analysed (scripts/30_donor_analyses.py). Also stores each unit's cell-class shares (all nuclei, unsorted).
Writes results/biology/donor_scores.tsv.gz and donor_composition.tsv.gz.
"""
from pathlib import Path
import importlib.util
import sys
import numpy as np
import pandas as pd
import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.contrasts import braak_value  # noqa: E402
from pdrep.harmonise import OUT as H  # noqa: E402

spec = importlib.util.spec_from_file_location("prog", ROOT / "scripts" / "27_programmes.py")
prog = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prog)
OUT = prog.OUT


def main():
    progs = prog.programmes()
    rows = []
    for ds_dir in sorted(prog.DE.iterdir()):
        units_f = ds_dir / "input" / "units.tsv"
        if not units_f.exists():
            continue
        units = pd.read_csv(units_f, sep="\t", dtype=str).set_index("unit")
        units["stage"] = braak_value(units.assign(braak_lb=pd.to_numeric(units.braak_lb, errors="coerce")))
        for f in sorted((ds_dir / "input").glob("*__*.csv.gz")):
            region, cls_f = f.name[: -len(".csv.gz")].split("__", 1)
            cls = cls_f.replace("_", " ").replace("Neurons other", "Neurons (other)")
            if cls not in prog.CLASSES:
                continue
            x = pd.read_csv(f, index_col=0)
            x.index = x.index.str.upper()
            lc = np.log2(x.div(x.sum(0), axis=1) * 1e6 + 1)
            z = lc.sub(lc.mean(1), axis=0).div(lc.std(1).replace(0, np.nan), axis=0)
            u = units.loc[z.columns]
            for pname, (genes, classes) in progs.items():
                if cls not in classes:
                    continue
                g = [x_ for x_ in genes if x_ in z.index]
                if len(g) < (2 if len(genes) == 2 else 5):
                    continue
                rows.append(pd.DataFrame({"dataset": ds_dir.name, "region": region, "unit": z.columns,
                                          "donor": u.donor.values, "disease": u.disease.values, "stage": u.stage.values,
                                          "age": pd.to_numeric(u.age, errors="coerce").values,
                                          "pmi": pd.to_numeric(u.pmi, errors="coerce").values, "sex": u.sex.values,
                                          "cell_class": cls, "programme": pname, "score": z.loc[g].mean(0).values}))
        print(ds_dir.name, flush=True)
    s = pd.concat(rows, ignore_index=True)
    s.to_csv(OUT / "donor_scores.tsv.gz", sep="\t", index=False)
    comp = []
    for f in sorted(H.glob("*.h5ad")):
        obs = ad.read_h5ad(f, backed="r").obs
        if "enrichment" in obs:
            obs = obs[~obs.enrichment.astype(str).str.lower().eq("sorted")]
        unit = obs.donor.astype(str) + "|" + obs.region.astype(str)
        c = pd.crosstab(unit, obs.cell_class.astype(str))
        share = c.div(c.sum(1), axis=0)
        neur = [k for k in c.columns if k in {"DA neurons", "MSNs", "Excitatory neurons", "Inhibitory neurons", "Neurons (other)"}]
        share["neurons_total"] = c[neur].sum(1) / c.sum(1) if neur else np.nan
        comp.append(share.reset_index().melt(id_vars="row_0" if "row_0" in share.reset_index() else share.index.name or "index",
                                             var_name="cell_class", value_name="share").assign(dataset=f.stem))
    cp = pd.concat(comp, ignore_index=True)
    cp = cp.rename(columns={cp.columns[0]: "unit"})
    cp.to_csv(OUT / "donor_composition.tsv.gz", sep="\t", index=False)
    print(s.groupby(["programme"]).size().to_string())


if __name__ == "__main__":
    main()

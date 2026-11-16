"""Which cell types change in proportion in PD? Random-effects meta-analysis across the harmonised datasets.

Usage: python -I scripts/22_composition_meta.py
Per dataset x region, donor units (sorted nuclei excluded) are compared case vs control on the logit of each cell
class's share, in two denominators: all nuclei, and non-neuronal nuclei (glial, vascular, immune classes; reduces
dependence on how much neuron-rich tissue a dissection captured). DA neurons are also expressed as a share of all
neurons. Effect = difference in mean logit share (case - control), SE from Welch. Regions of one dataset are combined
by fixed effect before pooling datasets with DerSimonian-Laird; region-group estimates pool dataset x region units.
Contrasts: PD/PDD vs control (main), any Lewy-body disease vs control (sensitivity).
Writes results/biology/composition_estimates.tsv and composition_meta.tsv.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
import anndata as ad

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.harmonise import OUT as H  # noqa: E402
from pdrep.contrasts import unit_groups  # noqa: E402
from pdrep.meta import fixed, dersimonian_laird, bh  # noqa: E402
from pdrep.plot.style import REGION_GROUP  # noqa: E402

from pdrep import config as cfg  # noqa: E402

OUT = cfg.RESULTS / "biology"
NEURONS = {"DA neurons", "MSNs", "Excitatory neurons", "Inhibitory neurons", "Neurons (other)"}
CONTRASTS = ["PDany_vs_control", "LBDany_vs_control"]


def welch(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    if len(a) < 3 or len(b) < 3:
        return np.nan, np.nan
    d = a.mean() - b.mean()
    se = np.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
    return d, se if se > 0 else np.nan


def logit_share(counts, total):
    p = (counts + 0.5) / (total + 1.0)
    return np.log(p / (1 - p))


def pool(est):
    """Per-dataset x region estimates -> pooled table: one estimate per dataset (regions combined by fixed effect),
    DerSimonian-Laird across datasets, per region group too; BH FDR per contrast x denominator x scope."""
    res = []
    for (con, cls, dname), d in est.groupby(["contrast", "cell_class", "denominator"]):
        # overall: one estimate per dataset (regions combined by fixed effect)
        per_ds = d.groupby("dataset").apply(lambda x: pd.Series(fixed(x.est, x.se), index=["est", "se"]))
        r = dersimonian_laird(per_ds.est, per_ds.se)
        res.append({"contrast": con, "cell_class": cls, "denominator": dname, "scope": "all regions", **r,
                    "datasets": ";".join(per_ds.dropna().index)})
        for rg, dd in d.groupby("region_group"):
            r = dersimonian_laird(dd.est, dd.se)
            res.append({"contrast": con, "cell_class": cls, "denominator": dname, "scope": rg, **r,
                        "datasets": ";".join(sorted(dd.dataset.unique()))})
    meta = pd.DataFrame(res)
    meta["fdr"] = np.nan
    for key, idx in meta.groupby(["contrast", "denominator", "scope"]).groups.items():
        meta.loc[idx, "fdr"] = bh(meta.loc[idx, "p"])
    return meta


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for f in sorted(H.glob("*.h5ad")):
        ds = f.stem
        obs = ad.read_h5ad(f, backed="r").obs
        if "enrichment" in obs:
            obs = obs[~obs.enrichment.astype(str).str.lower().eq("sorted")]
        obs = obs.assign(unit=obs.donor.astype(str) + "|" + obs.region.astype(str))
        counts = pd.crosstab(obs.unit, obs.cell_class.astype(str))
        meta = obs.groupby("unit", observed=True)[["donor", "region", "disease"] +
                                                   [c for c in ["braak_lb", "unified_lb_stage"] if c in obs]].first()
        g = unit_groups(meta)
        neur = [c for c in counts.columns if c in NEURONS]
        nonneur = [c for c in counts.columns if c not in NEURONS and c != "Unknown"]
        denoms = {"all nuclei": counts.sum(1), "non-neuronal nuclei": counts[nonneur].sum(1)}
        for reg in sorted(meta.region.astype(str).unique()):
            u = meta.index[meta.region.astype(str).eq(reg)]
            for con in CONTRASTS:
                case, ref = u[g.loc[u, con].eq("case")], u[g.loc[u, con].eq("ref")]
                if len(case) < 3 or len(ref) < 3:
                    continue
                for cls in counts.columns:
                    if cls == "Unknown":
                        continue
                    for dname, tot in denoms.items():
                        if dname == "non-neuronal nuclei" and cls in NEURONS:
                            continue
                        x = logit_share(counts[cls], tot)
                        d, se = welch(x[case], x[ref])
                        rows.append({"dataset": ds, "region": reg, "region_group": REGION_GROUP.get(reg, "other"),
                                     "contrast": con, "cell_class": cls, "denominator": dname, "est": d, "se": se,
                                     "n_case": len(case), "n_ref": len(ref),
                                     "share_case": float((counts.loc[case, cls] / tot[case]).mean()),
                                     "share_ref": float((counts.loc[ref, cls] / tot[ref]).mean())})
                    if cls == "DA neurons" and neur:
                        x = logit_share(counts[cls], counts[neur].sum(1))
                        d, se = welch(x[case], x[ref])
                        rows.append({"dataset": ds, "region": reg, "region_group": REGION_GROUP.get(reg, "other"),
                                     "contrast": con, "cell_class": cls, "denominator": "all neurons", "est": d, "se": se,
                                     "n_case": len(case), "n_ref": len(ref)})
        print(ds, "done", flush=True)
    est = pd.DataFrame(rows)
    est.to_csv(OUT / "composition_estimates.tsv", sep="\t", index=False)

    meta = pool(est)
    meta.to_csv(OUT / "composition_meta.tsv", sep="\t", index=False)
    show = meta[(meta.contrast == "PDany_vs_control") & (meta.scope == "all regions") & (meta.k >= 3)]
    print(show.sort_values(["denominator", "p"])[["denominator", "cell_class", "k", "est", "se", "p", "fdr", "I2"]]
          .round(3).to_string(index=False))


if __name__ == "__main__":
    main()

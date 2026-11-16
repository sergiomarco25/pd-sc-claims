"""Public variant of the patient axes (02_axes_interpretation.py, 03_axes_controls_clinical.py) without private data.

Usage: python -I scripts/subgroups/10_axes_public.py
The two consensus axes are rebuilt from the large cohorts that may be released (config [release] private_datasets
removed), on donor profiles without the gene modules of private papers ([release] private_modules). Every donor of
every public dataset is scored on them, and cases are compared with controls per dataset and pooled, exactly as in
03. Feeds the public release of F-axis2-shift ([release] rederived_findings, scripts/59b_public_findings.py).
Writes results/subgroups/public/axis_loadings.tsv, axis_scores.tsv, axes_shift.tsv and axes_vs_full.tsv (Tucker
congruence of each public axis with the full one, results/subgroups/axis_loadings.tsv).
"""
from pathlib import Path
import importlib.util
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pdrep import config as cfg  # noqa: E402

HERE = ROOT / "scripts" / "subgroups"
OUT = cfg.RESULTS / "subgroups" / "public"


def _module(name, file):
    spec = importlib.util.spec_from_file_location(name, HERE / file)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


ax2 = _module("axes02", "02_axes_interpretation.py")
ax3 = _module("axes03", "03_axes_controls_clinical.py")  # loads 01 as ax3.sg


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    sg = ax3.sg
    sg.STATES = [s for s in sg.STATES if s[0] not in cfg.PRIVATE_MODULES]
    cohorts = [c for c in ax2.COHORTS if cfg.is_public(c)]
    print("cohorts:", cohorts, "| modules left out:", sorted(cfg.PRIVATE_MODULES), flush=True)

    X, _ = sg.features()  # cases only (before all_features() adds the controls)
    load, _ = ax2.consensus_loadings(X, cohorts)
    load.to_csv(OUT / "axis_loadings.tsv", sep="\t")
    w = pd.DataFrame({"PC1": load[("mean", "PC1")], "PC2": load[("mean", "PC2")]})

    Xall, meta = ax3.all_features()
    s = ax3.axis_scores(Xall, w).merge(meta.reset_index()[["dataset", "donor", "disease"]], on=["dataset", "donor"])
    s = s[s.dataset.map(cfg.is_public)]
    s["group"] = np.where(s.disease.isin(ax3.CASES), "case", np.where(s.disease.eq("control"), "control", "other"))
    s.to_csv(OUT / "axis_scores.tsv", sep="\t", index=False)
    sh = ax3.shift_table(s)
    sh.to_csv(OUT / "axes_shift.tsv", sep="\t", index=False)
    print(sh.round(3).to_string(index=False), flush=True)

    full = cfg.RESULTS / "subgroups" / "axis_loadings.tsv"
    if full.exists():
        f = pd.read_csv(full, sep="\t", header=[0, 1], index_col=0)
        rows = []
        for axn in ("PC1", "PC2"):
            a, b = w[axn].align(f[("mean", axn)], join="inner")
            rows.append({"axis": axn, "features_public": len(w), "features_full": len(f), "features_shared": len(a),
                         "congruence": float(np.dot(a, b) / np.sqrt(np.dot(a, a) * np.dot(b, b)))})
        cg = pd.DataFrame(rows)
        cg.to_csv(OUT / "axes_vs_full.tsv", sep="\t", index=False)
        print(cg.round(3).to_string(index=False))


if __name__ == "__main__":
    main()

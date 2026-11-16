"""Rebuild results/de_cache/<dataset>/de_log.tsv from the written fit files (used when the R job died after the fits).

Usage: python -I scripts/12c_rebuild_de_log.py <dataset> [...]
n_case / n_ref come from input/manifest.tsv; braak fits get n_case = units with a stage value.
"""
from pathlib import Path
import sys
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pdrep.contrasts import braak_value  # noqa: E402

from pdrep import config as cfg  # noqa: E402

DE = cfg.RESULTS / "de_cache"
for ds in sys.argv[1:]:
    man = pd.read_csv(DE / ds / "input" / "manifest.tsv", sep="\t")
    units = pd.read_csv(DE / ds / "input" / "units.tsv", sep="\t").set_index("unit")
    rows = []
    for f in sorted((DE / ds / "de").glob("*.tsv.gz")):
        stem, con = f.name[:-len(".tsv.gz")].rsplit("__", 1)
        file = stem + ".csv.gz"
        t = pd.read_csv(f, sep="\t", usecols=["padj", "covariates"])
        m = man[(man.file == file) & (man.contrast == con)]
        reg, cls = (m.region.iloc[0], m.cell_class.iloc[0]) if len(m) else (stem.split("__")[0], man[man.file == file].cell_class.iloc[0])
        if con == "braak":
            cols = pd.read_csv(DE / ds / "input" / file, nrows=0).columns[1:]
            n_case, n_ref = int(braak_value(units.loc[cols]).notna().sum()), None
        else:
            n_case, n_ref = int(m.n_case.iloc[0]), int(m.n_ref.iloc[0])
        rows.append({"file": file, "region": reg, "cell_class": cls, "contrast": con, "n_case": n_case, "n_ref": n_ref,
                     "n_genes": len(t), "n_fdr05": int((t.padj < 0.05).sum()), "covariates": t.covariates.iloc[0]})
    pd.DataFrame(rows).to_csv(DE / ds / "de_log.tsv", sep="\t", index=False)
    print(ds, len(rows), "fits logged")

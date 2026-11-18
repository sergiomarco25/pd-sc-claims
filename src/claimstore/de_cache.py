"""DE cache of one dataset: edgeR pseudobulk fits for every region x cell class x contrast with enough units, and the
stage as a continuous covariate (src/claimstore/de_cache.R; inputs from claimstore.pseudobulk).

  python -I src/claimstore/cli.py de-cache DATASET [...]    run de_cache.R with this project's R ([envs] r_env) and
                                                            stage ([stage] in config/contrasts.toml)
  python -I src/claimstore/cli.py de-log DATASET [...]      rebuild results/de_cache/<dataset>/de_log.tsv from the
                                                            written fit files (when the R job died after the fits)

The stage contrast's code ([stage] contrast, default "stage") names its fit files (<region>__<class>__<code>.tsv.gz)
and its rows of de_log.tsv; claims use the same code for a continuous stage association.
Launching needs only the standard library; de-log needs pandas.
"""
import shutil
import subprocess
import sys
import tomllib
from pathlib import Path

from claimstore import config as cfg

HERE = Path(__file__).resolve().parent
DE = cfg.RESULTS / "de_cache"
CONTRASTS_FILE = cfg.REPO / "config" / "contrasts.toml"


def stage_spec():
    """[stage] of config/contrasts.toml (no pandas needed)."""
    if not CONTRASTS_FILE.is_file():
        return {}
    with open(CONTRASTS_FILE, "rb") as fh:
        return tomllib.load(fh).get("stage", {})


def stage_contrast(stage=None):
    st = stage_spec() if stage is None else stage
    return str(st.get("contrast", "stage")) if st.get("columns") else ""


def rscript():
    r = cfg.R_ENV / "bin" / "Rscript"
    return str(r) if r.is_file() else (shutil.which("Rscript") or "Rscript")


def command(dataset):
    """The Rscript command for one dataset."""
    st = stage_spec()
    cmd = [rscript(), "--vanilla", str(HERE / "de_cache.R"), str(DE / dataset)]
    code = stage_contrast(st)
    if code:
        cmd += [code, ",".join(st["columns"]), str(st.get("group_column", "")), str(st.get("zero_group", ""))]
    return cmd


def run(dataset):
    return subprocess.run(command(dataset), check=False).returncode


def rebuild_log(dataset):
    """de_log.tsv from the written fit files: n_case / n_ref from input/manifest.tsv; stage fits get n_case = units
    with a stage value."""
    import pandas as pd
    from claimstore.contrasts import stage_value
    code = stage_contrast()
    man = pd.read_csv(DE / dataset / "input" / "manifest.tsv", sep="\t")
    units = pd.read_csv(DE / dataset / "input" / "units.tsv", sep="\t").set_index("unit")
    rows = []
    for f in sorted((DE / dataset / "de").glob("*.tsv.gz")):
        stem, con = f.name[:-len(".tsv.gz")].rsplit("__", 1)
        file = stem + ".csv.gz"
        t = pd.read_csv(f, sep="\t", usecols=["padj", "covariates"])
        m = man[(man.file == file) & (man.contrast == con)]
        reg, cls = (m.region.iloc[0], m.cell_class.iloc[0]) if len(m) else (stem.split("__")[0], man[man.file == file].cell_class.iloc[0])
        if con == code:
            cols = pd.read_csv(DE / dataset / "input" / file, nrows=0).columns[1:]
            n_case, n_ref = int(stage_value(units.loc[cols]).notna().sum()), None
        else:
            n_case, n_ref = int(m.n_case.iloc[0]), int(m.n_ref.iloc[0])
        rows.append({"file": file, "region": reg, "cell_class": cls, "contrast": con, "n_case": n_case, "n_ref": n_ref,
                     "n_genes": len(t), "n_fdr05": int((t.padj < 0.05).sum()), "covariates": t.covariates.iloc[0]})
    pd.DataFrame(rows).to_csv(DE / dataset / "de_log.tsv", sep="\t", index=False)
    print(dataset, len(rows), "fits logged")


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0].startswith("-"):
        print("usage: de-cache DATASET [...]")
        return 2
    return max(run(ds) for ds in argv)


def rebuild_main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0].startswith("-"):
        print("usage: de-log DATASET [...]")
        return 2
    for ds in argv:
        rebuild_log(ds)
    return 0


if __name__ == "__main__":
    sys.exit(main())

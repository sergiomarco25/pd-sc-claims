"""Project paths and release settings, read from config/project.toml (claims/COMMUNITY_PLAN.md, step 1.1).

Every machine-specific path lives in one file, so the pipeline can run outside our cluster. In Python:

    from pdrep import config as cfg
    out = cfg.RESULTS / "biology"

Config file: $PDREP_CONFIG if set, else <repo>/config/project.toml. Each path can also be overridden by an
environment variable (PDREP_DATA_ROOT, PDREP_WORK, PDREP_LEGACY, PDREP_GENOME_FA, PDREP_PYTHON, PDREP_R_ENV). RAW, HARMONISED, RESULTS and
SCRATCH sit inside WORK. Values may use {repo} and {data_root}; relative paths are relative to the repository.

Shell and R scripts get the same paths as environment variables:

    eval "$(python3 -I src/pdrep/config.py --env)"      # exports PDREP_WORK, PDREP_RESULTS, ...
    python3 -I src/pdrep/config.py RESULTS              # prints one path

Standard library only (the MCP server container imports it through kg.py). Python 3.11+ reads the TOML file with
tomllib; older Pythons need `tomli`, otherwise only the environment variables and defaults apply.
"""
from __future__ import annotations

import os
import shlex
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]


def _read_toml(path: Path) -> dict:
    try:
        import tomllib as toml
    except ModuleNotFoundError:  # Python < 3.11
        try:
            import tomli as toml  # type: ignore[no-redef]
        except ModuleNotFoundError:
            print(f"pdrep.config: cannot read {path} (Python < 3.11 without tomli); using environment and defaults",
                  file=sys.stderr)
            return {}
    with open(path, "rb") as fh:
        return toml.load(fh)


def _load() -> tuple[dict, Path | None]:
    path = Path(os.environ.get("PDREP_CONFIG", REPO / "config" / "project.toml"))
    if path.is_file():
        return _read_toml(path), path
    return {}, None


CONFIG, CONFIG_FILE = _load()
_PATHS = CONFIG.get("paths", {})


def _path(key: str, default: str, **fmt) -> Path:
    raw = os.environ.get("PDREP_" + key.upper()) or _PATHS.get(key) or default
    raw = os.path.expanduser(os.path.expandvars(str(raw).format(repo=REPO, **fmt)))
    p = Path(raw)
    return p if p.is_absolute() else REPO / p


DATA_ROOT = _path("data_root", "data")                                    # project data area (large, not in git)
WORK = _path("work", "{data_root}", data_root=DATA_ROOT)                  # this pipeline: raw/, harmonised/, results/
LEGACY = _path("legacy", "{repo}/..")                                     # earlier project folders (revision/, metadata)
GENOME_FA = _path("genome_fa", "{data_root}/hg38.fa", data_root=DATA_ROOT)

# Analysis environments used by the SLURM runners (slurm/run_py.sbatch, slurm/run_r.sbatch).
_ENVS = CONFIG.get("envs", {})
PYTHON = Path(os.path.expanduser(os.environ.get("PDREP_PYTHON") or _ENVS.get("python") or sys.executable))
R_ENV = Path(os.path.expanduser(os.environ.get("PDREP_R_ENV") or _ENVS.get("r_env") or "."))

RAW = WORK / "raw"
HARMONISED = WORK / "harmonised"
RESULTS = WORK / "results"
SCRATCH = WORK / "scratch"
STORE_DB = RESULTS / "store" / "claim_store.sqlite"

# Pipeline settings (workflow/Snakefile): datasets, memory, contrasts.
PIPELINE = CONFIG.get("pipeline", {})
DATASETS = list(PIPELINE.get("datasets", []))


def mem_gb(dataset: str) -> int:
    """Memory for per-dataset jobs, from [pipeline] mem_gb / default_mem_gb."""
    return int(PIPELINE.get("mem_gb", {}).get(dataset, PIPELINE.get("default_mem_gb", 128)))


# Release rules (claims/COMMUNITY_PLAN.md, phase 3): datasets kept out of every public export.
RELEASE = CONFIG.get("release", {})
PRIVATE_DATASETS = frozenset(RELEASE.get("private_datasets", []))
PRIVATE_PAPERS = frozenset(RELEASE.get("private_papers", []))
SEPARABLE_FINDINGS = frozenset(RELEASE.get("separable_findings", []))
PRIVATE_MODULES = frozenset(RELEASE.get("private_modules", []))
REDERIVED_FINDINGS = dict(RELEASE.get("rederived_findings", {}))  # finding id -> {table, script}
PUBLIC_REPOSITORY = RELEASE.get("public_repository", "")


def is_public(dataset: str) -> bool:
    """True if results computed in this dataset may appear in a public export."""
    return dataset not in PRIVATE_DATASETS


PATHS = {"REPO": REPO, "DATA_ROOT": DATA_ROOT, "WORK": WORK, "LEGACY": LEGACY, "GENOME_FA": GENOME_FA,
         "RAW": RAW, "HARMONISED": HARMONISED, "RESULTS": RESULTS, "SCRATCH": SCRATCH, "STORE_DB": STORE_DB,
         "PYTHON": PYTHON, "R_ENV": R_ENV}


def _main(argv: list[str]) -> int:
    if not argv or argv[0] in ("-h", "--help"):
        print("usage: config.py --env | --show | NAME   (NAME in " + ", ".join(PATHS) + ")")
        return 0
    if argv[0] == "--env":
        for k, v in PATHS.items():
            print(f"export PDREP_{k}={shlex.quote(str(v))}")
        print(f"export PDREP_PRIVATE_DATASETS={shlex.quote(','.join(sorted(PRIVATE_DATASETS)))}")
        return 0
    if argv[0] == "--show":
        print(f"config file: {CONFIG_FILE or '(none: environment and defaults)'}")
        for k, v in PATHS.items():
            print(f"{k:11s} {v}")
        print(f"{'private':11s} {', '.join(sorted(PRIVATE_DATASETS)) or '-'}")
        return 0
    if argv[0] in PATHS:
        print(PATHS[argv[0]])
        return 0
    print(f"unknown name {argv[0]!r}", file=sys.stderr)
    return 2


if __name__ == "__main__":
    sys.exit(_main(sys.argv[1:]))

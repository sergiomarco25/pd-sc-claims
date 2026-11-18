"""Project configuration for the claim-store engine: paths, pipeline, release rules and project facts, read from one
TOML file (claims/COMMUNITY_PLAN.md steps 1.1 and 2.1).

Every machine-specific path and every fact about the project's field lives in the project's config file, so the
engine runs unchanged for another field. In Python:

    from claimstore import config as cfg
    out = cfg.RESULTS / "biology"

Config file: $CLAIMSTORE_CONFIG if set, else <repo>/config/project.toml. Each path can also be overridden by an
environment variable named <PREFIX>_<NAME> (DATA_ROOT, WORK, LEGACY, GENOME_FA, PYTHON, R_ENV, ...), where PREFIX is
`env_prefix` in the [project] table (default CLAIMSTORE). RAW, HARMONISED, RESULTS and SCRATCH sit inside WORK.
Values may use {repo} and {data_root}; relative paths are relative to the repository.

Shell and R scripts get the same paths as environment variables:

    eval "$(python3 -I src/claimstore/config.py --env)"      # exports <PREFIX>_WORK, <PREFIX>_RESULTS, ...
    python3 -I src/claimstore/config.py RESULTS              # prints one path

Project facts (tables of the config file):
  [project]  name, title, description, keywords, env_prefix, server_name: how releases, the server and the
             linked-data exports describe the store
  [claims]   paper_to_dataset: papers whose claims were made on another paper's cohort
  [rating]   enriched_sort (datasets whose composition is not tissue composition), shared_banks (datasets that may share
             donors), effect_templates / pattern_templates (test templates of each rating family), not_assessed (what
             the rating does not yet assess, stated on every rating)
  [server]   instructions (file with the server's instructions), examples (ids used in the tool descriptions)

Standard library only (the server container imports it through kg.py). Python 3.11+ reads the TOML file with tomllib;
older Pythons need `tomli`, otherwise only the environment variables and defaults apply.
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
            print(f"claimstore.config: cannot read {path} (Python < 3.11 without tomli); using environment and defaults",
                  file=sys.stderr)
            return {}
    with open(path, "rb") as fh:
        return toml.load(fh)


def _load() -> tuple[dict, Path | None]:
    path = Path(os.environ.get("CLAIMSTORE_CONFIG", REPO / "config" / "project.toml"))
    if path.is_file():
        return _read_toml(path), path
    return {}, None


CONFIG, CONFIG_FILE = _load()
PROJECT = CONFIG.get("project", {})
ENV_PREFIX = str(PROJECT.get("env_prefix", "CLAIMSTORE"))
_PATHS = CONFIG.get("paths", {})


def env(name: str, default: str | None = None) -> str | None:
    """The project's environment variable <PREFIX>_<NAME>."""
    return os.environ.get(f"{ENV_PREFIX}_{name}", default)


def _path(key: str, default: str, **fmt) -> Path:
    raw = env(key.upper()) or _PATHS.get(key) or default
    raw = os.path.expanduser(os.path.expandvars(str(raw).format(repo=REPO, **fmt)))
    p = Path(raw)
    return p if p.is_absolute() else REPO / p


DATA_ROOT = _path("data_root", "data")                                    # project data area (large, not in git)
WORK = _path("work", "{data_root}", data_root=DATA_ROOT)                  # this pipeline: raw/, harmonised/, results/
LEGACY = _path("legacy", "{repo}/..")                                     # earlier project folders read by a few scripts
GENOME_FA = _path("genome_fa", "{data_root}/hg38.fa", data_root=DATA_ROOT)

# Analysis environments used by the SLURM runners (slurm/run_py.sbatch, slurm/run_r.sbatch).
_ENVS = CONFIG.get("envs", {})
PYTHON = Path(os.path.expanduser(env("PYTHON") or _ENVS.get("python") or sys.executable))
R_ENV = Path(os.path.expanduser(env("R_ENV") or _ENVS.get("r_env") or "."))

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

# Facts about the project's claims, its rating and its server (see the module docstring).
CLAIMS = CONFIG.get("claims", {})
PAPER_TO_DATASET = dict(CLAIMS.get("paper_to_dataset", {}))
RATING = CONFIG.get("rating", {})
SERVER = CONFIG.get("server", {})


def is_public(dataset: str) -> bool:
    """True if results computed in this dataset may appear in a public export."""
    return dataset not in PRIVATE_DATASETS


def describe() -> dict:
    """How releases, the server and the linked-data exports name and describe the store ([project])."""
    return {"name": PROJECT.get("title", "Claim store"), "description": PROJECT.get("description", ""),
            "keywords": list(PROJECT.get("keywords", [])), "server_name": PROJECT.get("server_name", "claim-store")}


PATHS = {"REPO": REPO, "DATA_ROOT": DATA_ROOT, "WORK": WORK, "LEGACY": LEGACY, "GENOME_FA": GENOME_FA,
         "RAW": RAW, "HARMONISED": HARMONISED, "RESULTS": RESULTS, "SCRATCH": SCRATCH, "STORE_DB": STORE_DB,
         "PYTHON": PYTHON, "R_ENV": R_ENV}


def main(argv: list[str] | None = None) -> int:
    argv = sys.argv[1:] if argv is None else argv
    if not argv or argv[0] in ("-h", "--help"):
        print("usage: config.py --env | --show | NAME   (NAME in " + ", ".join(PATHS) + ")")
        return 0
    if argv[0] == "--env":
        for k, v in PATHS.items():
            print(f"export {ENV_PREFIX}_{k}={shlex.quote(str(v))}")
        print(f"export {ENV_PREFIX}_PRIVATE_DATASETS={shlex.quote(','.join(sorted(PRIVATE_DATASETS)))}")
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
    sys.exit(main())

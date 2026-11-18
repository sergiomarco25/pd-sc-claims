"""Static checks that keep the code portable (COMMUNITY_PLAN step 1.1). Standard library only.

- no hard-coded cluster paths in code (they belong in config/project.toml);
- every file that uses `cfg.` imports the config (pdrep.config or claimstore.config) before its first use (a multi-line statement once used it early).
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
def _is_file(p):
    try:
        return p.is_file()
    except OSError:  # unreadable (storage fault): reported by the other checks, not here
        return False


CODE = [p for d in ("scripts", "src", "workflow", "slurm", "deploy") for p in (ROOT / d).rglob("*")
        if p.suffix in {".py", ".sh", ".R", ".sbatch", ""} and "__pycache__" not in p.parts
        and "build" not in p.relative_to(ROOT).parts[:3] and _is_file(p)]
CLUSTER = re.compile(r"/(lustre|ictstr01)/")


def _text(p):
    try:
        return p.read_text()
    except (UnicodeDecodeError, OSError):
        return ""


def test_no_cluster_paths_in_code():
    hits = [str(p.relative_to(ROOT)) for p in CODE if p.name != "config.py" and CLUSTER.search(_text(p))]
    assert not hits, f"hard-coded cluster paths (use pdrep.config): {hits}"


def test_config_imported_before_use():
    bad = []
    for p in CODE:
        if p.suffix != ".py" or p.name == "config.py":
            continue
        lines = _text(p).split("\n")
        imp = next((i for i, ln in enumerate(lines)
                    if ln.startswith(("from pdrep import config as cfg", "from claimstore import config as cfg"))), None)
        use = next((i for i, ln in enumerate(lines) if "cfg." in ln and not ln.lstrip().startswith("#")), None)
        if use is not None and (imp is None or use < imp):
            bad.append(f"{p.relative_to(ROOT)}:{use + 1}")
    assert not bad, f"cfg used before `from pdrep import config as cfg`: {bad}"


if __name__ == "__main__":
    test_no_cluster_paths_in_code()
    test_config_imported_before_use()
    print("PASS")

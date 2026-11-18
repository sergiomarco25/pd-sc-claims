"""Rebuild results/de_cache/<dataset>/de_log.tsv from the written fit files (used when the R job died after the
fits): the engine's claimstore.de_cache.rebuild_log.

Run: python -I scripts/12c_rebuild_de_log.py <dataset> [...]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore.de_cache import *  # noqa: E402,F401,F403
from claimstore.de_cache import rebuild_main  # noqa: E402

if __name__ == "__main__":
    sys.exit(rebuild_main())

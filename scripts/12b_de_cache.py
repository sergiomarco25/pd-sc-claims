"""DE cache of one dataset (edgeR pseudobulk fits, src/claimstore/de_cache.R) with this project's R and Lewy-body
stage: the engine's claimstore.de_cache (see its docstring). Replaces scripts/12b_de_cache.R.

Run: python -I scripts/12b_de_cache.py <dataset> [...]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore.de_cache import *  # noqa: E402,F401,F403
from claimstore.de_cache import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())

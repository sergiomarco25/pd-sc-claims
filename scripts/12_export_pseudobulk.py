"""Export one dataset's pseudobulk for edgeR (counts per region x cell class, unit table with one group column per
contrast): the engine's claimstore.pseudobulk for this project (see its docstring).

Run: python -I scripts/12_export_pseudobulk.py <dataset>
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore.pseudobulk import *  # noqa: E402,F401,F403
from claimstore.pseudobulk import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())

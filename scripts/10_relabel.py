"""Re-run cell-class labelling, label check, donor table and pseudobulk for built datasets: the engine's
claimstore.harmonise.relabel for this project.

Run: python -I scripts/10_relabel.py <dataset> [<dataset> ...]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore.harmonise import *  # noqa: E402,F401,F403
from claimstore.harmonise import relabel_main  # noqa: E402

if __name__ == "__main__":
    sys.exit(relabel_main())

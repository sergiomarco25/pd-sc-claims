"""Re-run cell-class labelling, label check, donor table and pseudobulk for built datasets (pdrep.harmonise.relabel).

Usage: python -I scripts/10_relabel.py <dataset> [<dataset> ...]
"""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pdrep.harmonise import relabel  # noqa: E402

if __name__ == "__main__":
    for ds in sys.argv[1:]:
        relabel(ds)

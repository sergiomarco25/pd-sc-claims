"""Build one harmonised dataset: run its loader (src/pdrep/build/<dataset>.py), then pdrep.harmonise.finalise.

Usage: python -I scripts/09_build_dataset.py <dataset>
SLURM: sbatch -J h_<dataset> --mem=<GB>G slurm/run_py.sbatch scripts/09_build_dataset.py <dataset>
"""
from pathlib import Path
import importlib
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pdrep.harmonise import finalise  # noqa: E402

if __name__ == "__main__":
    ds = sys.argv[1]
    t0 = time.time()
    mod = importlib.import_module(f"pdrep.build.{ds}")
    a = mod.load()
    print(f"{ds}: loaded {a.shape} in {time.time() - t0:.0f}s", flush=True)
    finalise(a, ds, mod.AUTHOR_TO_CLASS)
    print(f"{ds}: done in {time.time() - t0:.0f}s")

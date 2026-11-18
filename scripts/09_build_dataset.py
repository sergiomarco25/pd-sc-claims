"""Build one harmonised dataset: run its loader (src/pdrep/build/<dataset>.py), then finalise: the engine's
claimstore.harmonise for this project (schema and labelling rules: config/harmonise.toml).

Run: python -I scripts/09_build_dataset.py <dataset>
SLURM: sbatch -J h_<dataset> --mem=<GB>G slurm/run_py.sbatch scripts/09_build_dataset.py <dataset>
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore.harmonise import *  # noqa: E402,F401,F403
from claimstore.harmonise import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())

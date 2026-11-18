"""Moved to the engine: claimstore.config (claims/COMMUNITY_PLAN.md step 2.1). This name stays for the PD scripts.

`import pdrep.config` and `from pdrep import config` give the claimstore module itself (same object, same state).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from claimstore import config as _engine  # noqa: E402

if __name__ == "__main__":  # `python3 -I src/pdrep/config.py --env | --show | NAME` (SLURM scripts, Snakefile)
    sys.exit(_engine.main())
sys.modules[__name__] = _engine

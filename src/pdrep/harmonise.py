"""Moved to the engine: claimstore.harmonise (claims/COMMUNITY_PLAN.md step 2.1); the PD schema, marker genes and
labelling rules are config/harmonise.toml. This name stays for the PD loaders and scripts.

`from pdrep.harmonise import RAW, OUT, finalise` gives the claimstore names themselves (same module object).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from claimstore import harmonise as _engine  # noqa: E402
sys.modules[__name__] = _engine

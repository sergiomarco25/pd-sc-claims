"""Claim region vocabulary -> harmonised dataset regions: the engine's (claimstore.contexts, config/contexts.toml)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from claimstore.contexts import REGION_MAP, target_regions  # noqa: E402,F401

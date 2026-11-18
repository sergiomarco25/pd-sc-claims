"""Moved to the engine: claimstore.store (claims/COMMUNITY_PLAN.md step 2.1). This name stays for the PD scripts.

`import pdrep.store` and `from pdrep import store` give the claimstore module itself (same object, same state).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from claimstore import store as _engine  # noqa: E402
sys.modules[__name__] = _engine

"""Moved to the engine: claimstore.h5io (claims/COMMUNITY_PLAN.md step 2.1). This name stays for the PD loaders.

`from pdrep.build._io import read_csr` gives the claimstore function itself (same module object).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from claimstore import h5io as _engine  # noqa: E402
sys.modules[__name__] = _engine

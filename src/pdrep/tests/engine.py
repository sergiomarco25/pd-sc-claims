"""Moved to the engine: claimstore.retest (claims/COMMUNITY_PLAN.md step 2.1); the PD GWAS test is pdrep.genetics,
called through [tests] genetic_test in config/project.toml. This name stays for the PD scripts.

`from pdrep.tests.engine import test_claim` gives the claimstore function itself (same module object, same state).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from claimstore import retest as _engine  # noqa: E402
sys.modules[__name__] = _engine

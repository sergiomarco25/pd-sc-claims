"""Moved to the engine: claimstore.mcp_server (claims/COMMUNITY_PLAN.md step 2.1). This name stays for the PD project.

`import pdrep.mcp_server` gives the claimstore module itself; running this file runs the engine's server
(`python3 -I src/pdrep/mcp_server.py`, as in .mcp.json and the tests).
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from claimstore import mcp_server as _engine  # noqa: E402

if __name__ == "__main__":
    _engine.main()
else:
    sys.modules[__name__] = _engine

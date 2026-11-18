"""Swap in the public tests of re-derived findings: the engine's claimstore.public_findings for this project (see its docstring).

Run: python3 -I scripts/59b_public_findings.py --db PUBLIC_STORE
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore.public_findings import *  # noqa: E402,F401,F403
from claimstore.public_findings import main  # noqa: E402

if __name__ == "__main__":
    main()

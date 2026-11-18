"""Compile the claim store (results/store/claim_store.sqlite and claims/store/): the engine's claimstore.build_store for this project (see its docstring).

Run: python -I scripts/51b_claim_store.py (analysis environment)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore.build_store import *  # noqa: E402,F401,F403
from claimstore.build_store import main  # noqa: E402

if __name__ == "__main__":
    main()

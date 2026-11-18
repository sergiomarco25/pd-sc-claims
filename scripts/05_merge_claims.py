"""Validate the per-paper claim files and merge them into claims/claims.yaml: the engine's claimstore.claims for this project (see its docstring).

Run: python -I scripts/05_merge_claims.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore.claims import *  # noqa: E402,F401,F403
from claimstore.claims import main  # noqa: E402

if __name__ == "__main__":
    main()

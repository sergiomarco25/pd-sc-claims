"""Cross-dataset meta-analysis of every claim: the engine's claimstore.claim_meta for this project (see its docstring
for the robustness classes). Other scripts load this file for `load`, `combine`, `robustness`, UP / DOWN and
RAW_SIGN_TEMPLATES.

Run: python -I scripts/15_meta_analysis.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore.claim_meta import *  # noqa: E402,F401,F403
from claimstore.claim_meta import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())

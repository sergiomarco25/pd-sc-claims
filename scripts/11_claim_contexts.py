"""Map every claim to its dataset, regions, cell classes and contrast (claims/contexts.tsv): the engine's claimstore.contexts for this project (see its docstring).

Run: python -I scripts/11_claim_contexts.py
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore.contexts import *  # noqa: E402,F401,F403
from claimstore.contexts import main  # noqa: E402

if __name__ == "__main__":
    main()

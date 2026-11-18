"""Public release of the claim store (claims/COMMUNITY_PLAN.md phase 4): the engine's claimstore.release for this project (see its docstring).

Run: python3 -I scripts/59_release.py VERSION [--previous OLD_STORE] [--out DIR] [--rating-python PATH] [--force]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore.release import *  # noqa: E402,F401,F403
from claimstore.release import main  # noqa: E402

if __name__ == "__main__":
    main()

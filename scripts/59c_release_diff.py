"""Per-claim changelog between two releases: the engine's claimstore.release_diff for this project (see its docstring).

Run: python3 -I scripts/59c_release_diff.py OLD_STORE NEW_STORE [--out DIR]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore.release_diff import *  # noqa: E402,F401,F403
from claimstore.release_diff import main  # noqa: E402

if __name__ == "__main__":
    main()

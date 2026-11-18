"""Linked-data exports of a release (JSON-LD, RO-Crate): the engine's claimstore.linked_data for this project (see its docstring).

Run: python3 -I scripts/59d_linked_data.py RELEASE_DIR
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore.linked_data import *  # noqa: E402,F401,F403
from claimstore.linked_data import main  # noqa: E402

if __name__ == "__main__":
    main()

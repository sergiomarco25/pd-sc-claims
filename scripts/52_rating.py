"""Evidence rating of every claim in the store (claims/RATING.md): the engine's claimstore.rating for this project.

Run: python -I scripts/52_rating.py [--db claim_store.sqlite] [--out folder]
     (defaults: the configured store and claims/store/; the release builder, scripts/59_release.py, re-rates a
     public store with --db/--out). Project facts used by the rating: [claims] and [rating] in config/project.toml.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore.rating import *  # noqa: E402,F401,F403
from claimstore.rating import main  # noqa: E402

if __name__ == "__main__":
    main()

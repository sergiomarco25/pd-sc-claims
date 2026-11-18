"""Audit a public atlas page for private material: the engine's claimstore.audit for this project (see its docstring).

Run: python3 -I scripts/61_audit_public.py [reports/pd_atlas_public.html]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore.audit import *  # noqa: E402,F401,F403
from claimstore.audit import main  # noqa: E402

if __name__ == "__main__" and len(sys.argv) < 2:  # this project's public atlas by default
    sys.argv.append(str(Path(__file__).resolve().parents[1] / "reports" / "pd_atlas_public.html"))

if __name__ == "__main__":
    main()

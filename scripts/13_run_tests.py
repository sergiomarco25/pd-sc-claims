"""Run claim tests in one dataset: the engine's claimstore.retest for this project (see its docstring; the PD GWAS
test is pdrep.genetics).

Run: python -I scripts/13_run_tests.py internal <dataset>    claims whose source dataset is <dataset>
     python -I scripts/13_run_tests.py cross <dataset>       all other papers' claims, in <dataset>
     add "gwas" as a third argument to run only GWAS-enrichment claims (written to <mode>_gwas/)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore.retest import *  # noqa: E402,F401,F403
from claimstore.retest import main  # noqa: E402

if __name__ == "__main__":
    sys.exit(main())

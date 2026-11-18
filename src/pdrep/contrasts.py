"""Contrasts of the PD project: the engine's (claimstore.contrasts, defined in config/contrasts.toml) under PD names.

LEWY is the `lewy_body` group of the config file; braak_value is the engine's stage (Braak Lewy-body stage, else the
Unified LB stage, with pathology-free controls at 0). `braak` (the stage as a continuous covariate) is not a group
contrast and is handled by the DE cache and the tests engine directly.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from claimstore import contrasts as _engine  # noqa: E402
from claimstore.contrasts import CONTRASTS, GROUPS, SPEC_FILE, STAGE, unit_groups  # noqa: E402,F401

LEWY = GROUPS.get("lewy_body", set())


def braak_value(meta):
    """Lewy-body stage per unit: braak_lb, else the Unified LB stage, else NaN; controls without a stage count as
    stage 0 when the dataset stages its cases."""
    return _engine.stage_value(meta)


def load(path=SPEC_FILE):
    """config/contrasts.toml -> (CONTRASTS, LEWY), as before the engine split."""
    contrasts, groups, _ = _engine.load(path, default_stage=STAGE)
    return contrasts, groups.get("lewy_body", set())

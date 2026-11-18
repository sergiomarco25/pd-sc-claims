"""Contrasts: which donor x region units are case / reference, defined as data in config/contrasts.toml.

Every group contrast is defined once (format in that file), so every analysis uses the same groups and a new contrast
needs no code. A side is a table of conditions on unit metadata (equal, one of, contains, stage range); the stage
itself ([stage] in the file) is the first of its columns with values, with the reference group counted as stage 0.
Top-level lists of the file are named groups (GROUPS).
"""
import tomllib

import numpy as np
import pandas as pd

from claimstore import config as cfg

SPEC_FILE = cfg.REPO / "config" / "contrasts.toml"
_RANGE_KEYS = {"min", "max"}


def _col(m, c):
    return m[c].astype(str) if c in m else pd.Series("", index=m.index)


def stage_value(meta, stage=None):
    """Stage per unit: the first configured column with any value; units of the zero group without a stage count as
    stage 0 when anyone in the table has one."""
    st = stage if stage is not None else STAGE
    b = pd.Series(np.nan, index=meta.index)
    for c in st.get("columns", []):
        if c in meta:
            v = pd.to_numeric(meta[c], errors="coerce")
            if v.notna().any():
                b = v
                break
    if b.notna().any() and st.get("zero_group") is not None:
        grp = _col(meta, st.get("group_column", "group"))
        b = b.where(~(b.isna() & grp.eq(str(st["zero_group"]))), 0.0)
    return b


def _condition(column, value, where, stage):
    """One condition of a contrast side -> function(meta) -> boolean Series."""
    if column == stage.get("name"):
        if not isinstance(value, dict) or not value or set(value) - _RANGE_KEYS:
            raise ValueError(f"{where}: {column} needs {{min = .., max = ..}}, got {value!r}")
        lo, hi = value.get("min", -np.inf), value.get("max", np.inf)
        return lambda m: stage_value(m, stage).between(lo, hi)
    if isinstance(value, dict):
        if set(value) != {"contains"}:
            raise ValueError(f"{where}: {column} = {value!r}; a table condition must be {{contains = \"...\"}}")
        x = str(value["contains"]).lower()
        return lambda m: _col(m, column).str.lower().str.contains(x, regex=False)
    if isinstance(value, list):
        vals = [str(v) for v in value]
        return lambda m: _col(m, column).isin(vals)
    if isinstance(value, (str, int, float)) and not isinstance(value, bool):
        v = str(value)
        return lambda m: _col(m, column).eq(v)
    raise ValueError(f"{where}: {column} = {value!r} is not a text, a list or {{contains = ...}}")


def _side(conditions, where, stage):
    if not isinstance(conditions, dict) or not conditions:
        raise ValueError(f"{where}: needs at least one condition")
    parts = [_condition(c, v, where, stage) for c, v in conditions.items()]

    def select(m):
        out = pd.Series(True, index=m.index)
        for p in parts:
            out &= p(m).fillna(False).astype(bool)
        return out
    return select


def load(path=SPEC_FILE, default_stage=None):
    """contrasts.toml -> (CONTRASTS {name: (case selector, reference selector)}, GROUPS {name: set}, STAGE dict).
    A file without a [stage] table uses `default_stage`."""
    with open(path, "rb") as fh:
        spec = tomllib.load(fh)
    stage = spec.get("stage", default_stage or {})
    contrasts = {}
    for name, c in spec.get("contrasts", {}).items():
        if set(c) != {"case", "ref"}:
            raise ValueError(f"contrast {name}: needs exactly `case` and `ref`, got {sorted(c)}")
        contrasts[name] = (_side(c["case"], f"{name}.case", stage), _side(c["ref"], f"{name}.ref", stage))
    groups = {k: set(v) for k, v in spec.items() if isinstance(v, list)}
    return contrasts, groups, stage


CONTRASTS, GROUPS, STAGE = load()


def unit_groups(meta, group_column=None):
    """DataFrame (units x contrasts) with 'case', 'ref' or '' per unit."""
    m = meta.copy()
    gc = group_column or STAGE.get("group_column", "group")
    if gc in m:
        m[gc] = m[gc].astype(str)
    out = {}
    for name, (case, ref) in CONTRASTS.items():
        c, r = case(m).fillna(False).values, ref(m).fillna(False).values
        out[name] = np.where(c & ~r, "case", np.where(r & ~c, "ref", ""))
    return pd.DataFrame(out, index=meta.index)

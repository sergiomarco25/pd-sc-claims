"""Contrast definitions: which donor x region units are case / reference for each contrast.

Every contrast in scripts/11_claim_contexts.py that compares groups of donors is defined once, as data in
config/contrasts.toml (format described there), so the DE cache, the composition tests and the cross-dataset tests
all use the same groups and a new contrast needs no code. `braak` is a continuous covariate (Lewy-body Braak stage,
or the Unified LB stage where Braak is absent) and is handled separately.
"""
import tomllib
from pathlib import Path

import numpy as np
import pandas as pd

SPEC_FILE = Path(__file__).resolve().parents[2] / "config" / "contrasts.toml"
_RANGE_KEYS = {"min", "max"}


def _col(m, c):
    return m[c].astype(str) if c in m else pd.Series("", index=m.index)


def braak_value(meta):
    """Lewy-body stage per unit: braak_lb, else the Unified LB stage (menon2026), else NaN.

    Controls without a stage count as stage 0 when the dataset stages its cases (pathology-free by definition)."""
    b = pd.to_numeric(meta.get("braak_lb"), errors="coerce") if "braak_lb" in meta else pd.Series(np.nan, index=meta.index)
    if b.isna().all() and "unified_lb_stage" in meta:
        b = pd.to_numeric(meta.unified_lb_stage, errors="coerce")
    if b.notna().any():
        b = b.where(~(b.isna() & meta.disease.astype(str).eq("control")), 0.0)
    return b


def _condition(column, value, where):
    """One condition of a contrast side -> function(meta) -> boolean Series."""
    if column == "lewy_stage":
        if not isinstance(value, dict) or not value or set(value) - _RANGE_KEYS:
            raise ValueError(f"{where}: lewy_stage needs {{min = .., max = ..}}, got {value!r}")
        lo, hi = value.get("min", -np.inf), value.get("max", np.inf)
        return lambda m: braak_value(m).between(lo, hi)
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


def _side(conditions, where):
    if not isinstance(conditions, dict) or not conditions:
        raise ValueError(f"{where}: needs at least one condition")
    parts = [_condition(c, v, where) for c, v in conditions.items()]

    def select(m):
        out = pd.Series(True, index=m.index)
        for p in parts:
            out &= p(m).fillna(False).astype(bool)
        return out
    return select


def load(path=SPEC_FILE):
    """config/contrasts.toml -> (CONTRASTS {name: (case selector, reference selector)}, LEWY set)."""
    with open(path, "rb") as fh:
        spec = tomllib.load(fh)
    contrasts = {}
    for name, c in spec.get("contrasts", {}).items():
        if set(c) != {"case", "ref"}:
            raise ValueError(f"contrast {name}: needs exactly `case` and `ref`, got {sorted(c)}")
        contrasts[name] = (_side(c["case"], f"{name}.case"), _side(c["ref"], f"{name}.ref"))
    return contrasts, set(spec.get("lewy_body", []))


CONTRASTS, LEWY = load()


def unit_groups(meta):
    """DataFrame (units x contrasts) with 'case', 'ref' or '' per unit."""
    m = meta.copy()
    m["disease"] = m.disease.astype(str)
    out = {}
    for name, (case, ref) in CONTRASTS.items():
        c, r = case(m).fillna(False).values, ref(m).fillna(False).values
        out[name] = np.where(c & ~r, "case", np.where(r & ~c, "ref", ""))
    return pd.DataFrame(out, index=meta.index)

"""Swap in the public tests of re-derived findings (claims/COMMUNITY_PLAN.md, phase 4).

Usage: python -I scripts/59b_public_findings.py --db PUBLIC_STORE
Called by scripts/59_release.py after private material is removed and before re-rating. For every finding in config
[release] rederived_findings, its evidence table is replaced by the public one (computed without private data, e.g.
by scripts/subgroups/10_axes_public.py), private datasets are dropped from its filter, and its per-dataset tests
(claim_dataset) are recomputed with claimstore.store.finding_tests, as the store build does. Its claim row then
points at the public script and table. Needs pandas, scipy and pyyaml. Prints a JSON summary.
"""
import argparse
import copy
import json
import math
import sqlite3
import sys
from pathlib import Path

from claimstore import config as cfg
from claimstore import store

ROOT = cfg.REPO


def _sql(v):
    if v is None:
        return None
    if hasattr(v, "item"):  # numpy scalar
        v = v.item()
    return None if isinstance(v, float) and math.isnan(v) else v


def rederive(con, findings, rederived, is_public=cfg.is_public):
    """Replace the claim_dataset rows of each re-derived finding by tests from its public table.
    findings: {id: findings.yaml entry}; rederived: {id: {"table": path under results/, "script": path}}."""
    cols = [r[1] for r in con.execute("PRAGMA table_info(claim_dataset)")]
    ccols = {r[1] for r in con.execute("PRAGMA table_info(claim)")}
    summary = {}
    for fid, spec in sorted(rederived.items()):
        if fid not in findings:
            raise KeyError(f"{fid}: not in claims/findings.yaml")
        f = copy.deepcopy(findings[fid])
        ev = f["evidence"]
        ev["table"] = spec["table"]
        flt = ev.get("filter") or {}
        if isinstance(flt.get("dataset"), list):
            flt["dataset"] = [d for d in flt["dataset"] if is_public(d)]
        res = store.resolve_evidence(ev)
        res = res[res.dataset.map(is_public)] if len(res) else res
        if res.empty:
            raise ValueError(f"{fid}: no public rows in results/{spec['table']}")
        tests = store.finding_tests(f, res)
        con.execute("DELETE FROM claim_dataset WHERE claim_id = ?", (fid,))
        con.executemany(f"INSERT INTO claim_dataset ({', '.join(cols)}) VALUES ({', '.join('?' * len(cols))})",
                        [[_sql(t.get(c)) for c in cols] for t in tests])
        upd = {"analysis_script": spec.get("script", ""), "analysis_results": f"results/{spec['table']}",
               "location": f"{spec.get('script', '')} -> results/{spec['table']}"}
        upd = {k: v for k, v in upd.items() if k in ccols}
        if upd:
            con.execute(f"UPDATE claim SET {', '.join(k + ' = ?' for k in upd)} WHERE claim_id = ?", [*upd.values(), fid])
        summary[fid] = {"table": spec["table"], "script": spec.get("script", ""), "tests": len(tests),
                        "datasets": sorted(str(d) for d in res.dataset)}
    con.commit()
    return summary


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--db", required=True, help="public store being built (modified in place)")
    a = ap.parse_args()
    import yaml
    findings = {f["id"]: f for f in yaml.safe_load(store.FINDINGS.read_text()) or []}
    con = sqlite3.connect(a.db)
    present = {r[0] for r in con.execute("SELECT claim_id FROM claim WHERE claim_id LIKE 'F-%'")}
    todo = {k: v for k, v in cfg.REDERIVED_FINDINGS.items() if k in present}
    print(json.dumps(rederive(con, findings, todo)))
    con.close()


if __name__ == "__main__":
    main()

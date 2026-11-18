"""Per-claim changelog between two claim-store releases (claims/COMMUNITY_PLAN.md, phase 4 step 5).

Usage: python3 -I scripts/59c_release_diff.py OLD_STORE NEW_STORE [--out DIR]
Compares two stores claim by claim (e.g. results/releases/v1.0.0/claim_store.sqlite and v1.1.0/) and writes
DIR/changes.tsv, one row per change (claim_id, change, old, new, statement), and DIR/CHANGES.md: counts, then every
change that alters what a reader would conclude (verdict, certainty, claims added or removed). Changes recorded:
  added, removed           claim present in one store only
  verdict, certainty       rating changed
  cohorts                  number of independent cohorts in the rating
  effect                   pooled effect moved by 0.05 or more (approx. SMD)
  statement                claim text edited
  datasets_added/removed   datasets in which the claim is tested
  dataset_added/removed    datasets of the store (claim_id empty)
Standard library only; scripts/59_release.py runs it with --previous.
"""
import argparse
import csv
import sqlite3
from pathlib import Path

EFFECT_STEP = 0.05
FIELDS = ["claim_id", "change", "old", "new", "statement"]


def _claims(con):
    cols = {r[1] for r in con.execute("PRAGMA table_info(claim)")}
    want = [c for c in ("claim_id", "statement", "verdict", "certainty", "k_cohorts", "est") if c in cols]
    rows = {r[0]: dict(zip(want, r)) for r in con.execute(f"SELECT {', '.join(want)} FROM claim")}
    tested = {}
    for cid, ds in con.execute("SELECT claim_id, dataset FROM claim_dataset WHERE role IS NULL OR role != 'source'"):
        tested.setdefault(cid, set()).add(ds)
    for cid, r in rows.items():
        r["datasets"] = tested.get(cid, set())
    return rows


def label(con, path):
    try:
        return con.execute("SELECT [commit] FROM build_info").fetchone()[0] or str(path)
    except sqlite3.Error:
        return str(path)


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def diff(old_db, new_db):
    """Return (changes as dicts, labels (old, new))."""
    a = sqlite3.connect(f"file:{old_db}?mode=ro", uri=True)
    b = sqlite3.connect(f"file:{new_db}?mode=ro", uri=True)
    old, new = _claims(a), _claims(b)
    rows = []

    def add(cid, change, o, n, statement=""):
        rows.append({"claim_id": cid, "change": change, "old": "" if o is None else o, "new": "" if n is None else n,
                     "statement": statement or ""})

    for cid in sorted(set(old) | set(new)):
        o, n = old.get(cid), new.get(cid)
        if o is None:
            add(cid, "added", "", n.get("verdict"), n.get("statement"))
            continue
        if n is None:
            add(cid, "removed", o.get("verdict"), "", o.get("statement"))
            continue
        st = n.get("statement")
        for f in ("verdict", "certainty"):
            if (o.get(f) or "") != (n.get(f) or ""):
                add(cid, f, o.get(f), n.get(f), st)
        if _num(o.get("k_cohorts")) != _num(n.get("k_cohorts")):
            add(cid, "cohorts", o.get("k_cohorts"), n.get("k_cohorts"), st)
        eo, en = _num(o.get("est")), _num(n.get("est"))
        if (eo is None) != (en is None) or (eo is not None and abs(eo - en) >= EFFECT_STEP):
            add(cid, "effect", "" if eo is None else f"{eo:.2f}", "" if en is None else f"{en:.2f}", st)
        if (o.get("statement") or "") != (st or ""):
            add(cid, "statement", o.get("statement"), st, st)
        if n["datasets"] - o["datasets"]:
            add(cid, "datasets_added", "", ";".join(sorted(n["datasets"] - o["datasets"])), st)
        if o["datasets"] - n["datasets"]:
            add(cid, "datasets_removed", ";".join(sorted(o["datasets"] - n["datasets"])), "", st)
    dso = {r[0] for r in a.execute("SELECT dataset FROM dataset")}
    dsn = {r[0] for r in b.execute("SELECT dataset FROM dataset")}
    for ds in sorted(dsn - dso):
        add("", "dataset_added", "", ds)
    for ds in sorted(dso - dsn):
        add("", "dataset_removed", ds, "")
    labels = (label(a, old_db), label(b, new_db))
    a.close()
    b.close()
    return rows, labels


def write(rows, labels, out):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    with open(out / "changes.tsv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    counts = {}
    for r in rows:
        counts[r["change"]] = counts.get(r["change"], 0) + 1
    md = [f"# Changes from {labels[0]} to {labels[1]}", "",
          "Every change, one row per claim and kind, is in `changes.tsv`.", "", "| change | claims |", "|---|---|"]
    md += [f"| {k} | {v} |" for k, v in sorted(counts.items())] or ["| (none) | 0 |"]
    for kind, title in (("verdict", "Verdict changed"), ("certainty", "Certainty changed")):
        sel = [r for r in rows if r["change"] == kind]
        if sel:
            md += ["", f"## {title}", ""]
            md += [f"- `{r['claim_id']}`: {r['old'] or '(none)'} → {r['new'] or '(none)'}. {_short(r['statement'])}" for r in sel]
    for kind, title in (("added", "Claims added"), ("removed", "Claims removed")):
        sel = [r["claim_id"] for r in rows if r["change"] == kind]
        if sel:
            by = {}
            for cid in sel:
                by.setdefault(cid.rsplit("-", 1)[0] if not cid.startswith(("F-", "BG-")) else cid.split("-")[0], []).append(cid)
            md += ["", f"## {title} ({len(sel)})", ""]
            md += [f"- {k}: {', '.join(f'`{c}`' for c in v)}" for k, v in sorted(by.items())]
    (out / "CHANGES.md").write_text("\n".join(md) + "\n")
    return counts


def _short(s, n=140):
    s = (s or "").strip()
    return s if len(s) <= n else s[: n - 1].rstrip() + "…"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("old")
    ap.add_argument("new")
    ap.add_argument("--out", default=".", help="folder for changes.tsv and CHANGES.md")
    a = ap.parse_args()
    rows, labels = diff(a.old, a.new)
    counts = write(rows, labels, a.out)
    print(f"{labels[0]} -> {labels[1]}: {counts or 'no changes'}")


if __name__ == "__main__":
    main()

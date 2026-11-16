"""Build the small test store tests/fixture/claim_store.sqlite from the full claim store (COMMUNITY_PLAN step 1.5).

The fixture is committed and used by CI, so it must hold nothing private: claims of private datasets' papers, results
computed in private datasets (config [release] private_datasets) and this study's findings are left out, questions
are emptied and consensus statements keep their wording but not our evidence notes. Aggregate ratings of the kept
claims still reflect every dataset; the public release (COMMUNITY_PLAN phase 4) recomputes them without private ones.

Usage: python3 -I tests/fixture/make_fixture.py [--db full_store.sqlite] [--per-paper 3]
Standard library only.
"""
import argparse
import sqlite3
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from pdrep import config as cfg  # noqa: E402

OUT = Path(__file__).resolve().parent / "claim_store.sqlite"
ALWAYS = ["xu2023-C11"]  # used by tests/test_mcp_http.py
MAX_GENES = 200          # per signature


def cols(con, table):
    return [r[1] for r in con.execute(f"PRAGMA table_info('{table}')")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", default=str(cfg.STORE_DB))
    ap.add_argument("--per-paper", type=int, default=3)
    a = ap.parse_args()

    private = sorted(cfg.PRIVATE_DATASETS)
    src = sqlite3.connect(f"file:{a.db}?mode=ro", uri=True)
    OUT.unlink(missing_ok=True)
    con = sqlite3.connect(OUT)
    src.backup(con)
    src.close()

    ph = ",".join("?" * len(private)) or "''"
    private_papers = [r[0] for r in con.execute(f"SELECT DISTINCT paper FROM dataset WHERE dataset IN ({ph})", private)]
    papers = {p for d in private_papers for p in d.split(";")}
    keep = set(ALWAYS)
    rows = con.execute("SELECT claim_id, paper, source_dataset, cell_types FROM claim ORDER BY claim_id").fetchall()
    per = {}
    for cid, paper, sds, cts in rows:
        if cid.startswith("F-") or paper in papers or (sds or "") in private:
            continue
        if per.get(paper, 0) < a.per_paper or ("Microglia" in (cts or "") and per.get(("mg", paper), 0) < 1):
            keep.add(cid)
            per[paper] = per.get(paper, 0) + 1
            if "Microglia" in (cts or ""):
                per[("mg", paper)] = 1
    con.execute("CREATE TEMP TABLE keep(claim_id TEXT PRIMARY KEY)")
    con.executemany("INSERT INTO keep VALUES (?)", [(c,) for c in sorted(keep)])

    tables = [r[0] for r in con.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")]
    for t in tables:
        c = cols(con, t)
        if "claim_id" in c:
            con.execute(f"DELETE FROM '{t}' WHERE claim_id NOT IN (SELECT claim_id FROM keep)")
        if "dataset" in c and private:
            con.execute(f"DELETE FROM '{t}' WHERE dataset IN ({ph})", private)
    con.execute("DELETE FROM claim_relation WHERE claim_a NOT IN (SELECT claim_id FROM keep) "
                "OR claim_b NOT IN (SELECT claim_id FROM keep)")
    con.execute("DELETE FROM signature WHERE signature_id NOT IN (SELECT signature_id FROM claim_signature)")
    con.execute("DELETE FROM signature_gene WHERE signature_id NOT IN (SELECT signature_id FROM signature) "
                f"OR CAST(rank AS INTEGER) > {MAX_GENES}")
    con.execute("UPDATE consensus_statement SET our_evidence = ''")
    con.execute("DELETE FROM question")
    con.execute("UPDATE build_info SET [commit] = 'test fixture', n_claims = (SELECT count(*) FROM claim), "
                "n_datasets = (SELECT count(*) FROM dataset)")
    con.commit()
    con.execute("VACUUM")
    n = {t: con.execute(f"SELECT count(*) FROM '{t}'").fetchone()[0] for t in ("claim", "claim_dataset", "dataset")}
    con.close()
    print(f"{OUT} ({OUT.stat().st_size / 1e6:.1f} MB): {n['claim']} claims, {n['claim_dataset']} tests, "
          f"{n['dataset']} datasets; private datasets removed: {', '.join(private) or '-'}")


if __name__ == "__main__":
    main()

"""MCP server for a claim store and its knowledge graph (claims/STORE_DESIGN.md section 3, "AI access").

Lets Claude (Claude Code, Claude Desktop, claude.ai via a connector) or any other MCP client answer from the store and
cite claim ids instead of answering from memory.

Transport: stdio (newline-delimited JSON-RPC 2.0), implemented with the standard library only, so it runs with any
Python 3.9+ and needs no install (for the remote, Streamable HTTP + REST/OpenAPI version see mcp_http.py):
    python3 -I src/claimstore/mcp_server.py      (or the project's entry point)
From another machine, run it over ssh:
    {"command": "ssh", "args": ["<host>", "python3 -I <repo>/src/claimstore/mcp_server.py"]}
Name, title, instructions and the example ids in the tool descriptions come from the project config ([project],
[server]; claimstore.config).

Environment (<PREFIX> = [project] env_prefix):
  <PREFIX>_STORE_DB        path of claim_store.sqlite (default: the project results/store copy)
  <PREFIX>_MCP_READONLY=1  disable add_question (the only tool that writes: it appends to claims/questions.yaml)

The store is opened read-only; the graph is rebuilt when the database file changes (after 51b_claim_store.py).
"""
import json
import os
import re
import sqlite3
import sys
import threading
import time
import traceback
from datetime import date
from pathlib import Path

if __package__ in (None, ""):  # run as a file: make the package importable (repository src/, or /app in a container)
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from claimstore import config as cfg  # noqa: E402
from claimstore import kg  # noqa: E402

ROOT = cfg.REPO
DB = Path(cfg.env("STORE_DB") or kg.DEFAULT_DB)
READONLY = cfg.env("MCP_READONLY") == "1"
QUESTIONS = ROOT / "claims" / "questions.yaml"
PROTOCOLS = ["2025-11-25", "2025-06-18", "2025-03-26", "2024-11-05"]
VERSION = "0.2.0"
SQL_SECONDS = float(cfg.env("SQL_SECONDS", "5"))  # wall-clock limit per sql query
MAX_ROWS = 1000  # hard cap on any `limit` argument

GENERIC_INSTRUCTIONS = """\
Claim store: published conclusions of a field's studies, extracted as structured claims, each re-tested in its source \
dataset and in every other harmonised dataset, then rated (GRADE-like: verdict + certainty). Also holds findings of \
the study that built it (origin=finding, ids F-...), consensus statements (CS..), reference knowledge (BG-...) and \
questions.

How to use:
- Start with `overview`; find records with `search` or `search_claims`; read one with `get_claim`, `get_statement`, \
`get_dataset` or `gene`; explore links with `neighbours` and `connect` (knowledge graph); anything else with `sql`.
- Cite claim ids and statement ids for every statement you make from the store.
- Report the certainty with the verdict; 'Untested' means no dataset could test it, not that it is false.
- If the store cannot answer a question, say so; `add_question` logs it for a new analysis.
"""


def _instructions():
    """The project's server instructions ([server] instructions: a file next to the config file), else generic."""
    f = cfg.SERVER.get("instructions")
    if f and cfg.CONFIG_FILE and (cfg.CONFIG_FILE.parent / f).is_file():
        return (cfg.CONFIG_FILE.parent / f).read_text()
    return GENERIC_INSTRUCTIONS


INSTRUCTIONS = _instructions()
ABOUT = cfg.describe()
# ids and names quoted in the tool descriptions ([server] examples)
EX = {"claim": "paper2024-C01", "finding": "F-example", "statement": "CS01", "gene": "GENE1", "cell": "CellType",
      "region": "REGION", "dataset": "dataset2024", "background": "BG-example", "paper": "paper2024",
      "cell_types": "a cell type, All", "regions": "a region", **cfg.SERVER.get("examples", {})}


# ------------------------------------------------------------------ store access
class Store:
    def __init__(self, path):
        self.path = path
        self._mtime = None
        self._g = None
        self._lock = threading.Lock()

    def con(self, seconds=None):
        c = kg.connect(self.path)
        c.row_factory = sqlite3.Row
        c.execute("PRAGMA query_only = 1")
        if seconds:  # abort runaway queries (recursive CTEs, cross joins) from public clients
            end = time.monotonic() + seconds
            c.set_progress_handler(lambda: int(time.monotonic() > end), 10000)
        return c

    def rows(self, sql, args=()):
        with self.con() as c:
            return [dict(r) for r in c.execute(sql, args)]

    @property
    def g(self):
        m = self.path.stat().st_mtime
        if self._g is None or m != self._mtime:
            with self._lock:
                if self._g is None or m != self._mtime:
                    self._g, self._mtime = kg.build(self.path), m
        return self._g


S = Store(DB)


def _clean(d, drop=()):
    return {k: v for k, v in d.items() if v not in (None, "", "[]", "{}") and k not in drop}


def _node_brief(nid):
    n = S.g.nodes[nid]
    out = dict(id=nid, type=n["type"], label=n["label"])
    for k in ("verdict", "certainty", "paper", "origin", "short", "n_claims", "n_donors", "status"):
        if k in n:
            out[k] = n[k]
    return out


def _resolve_one(name, types=None):
    hits = S.g.resolve(name, types)
    if not hits:
        raise ToolError(f"No node matches {name!r}. Use `search` to find ids (e.g. claim:{EX['claim']}, gene:{EX['gene']}, "
                        f"cell:{EX['cell']}, region:{EX['region']}, dataset:{EX['dataset']}, statement:{EX['statement']}).")
    if len(hits) > 1 and not (S.g.nodes[hits[0]]["label"].lower() == str(name).lower() or hits[0].split(":", 1)[1]
                              .lower() == str(name).lower() or hits[0] == name):
        raise ToolError(f"{name!r} is ambiguous: " + ", ".join(hits[:15]))
    return hits[0]


CLAIM_BRIEF = ("claim_id, paper, origin, statement, cell_types, regions, contrast, direction, claim_type, verdict, "
               "certainty, k_cohorts, est, lo, hi")


# ------------------------------------------------------------------ tools
class ToolError(Exception):
    pass


def _int(v, default, lo, hi, name="limit"):
    try:
        v = default if v in (None, "") else int(v)
    except (TypeError, ValueError):
        raise ToolError(f"{name} must be an integer") from None
    return max(lo, min(v, hi))


def t_overview():
    g = S.g
    info = S.rows("SELECT * FROM build_info")
    verdicts = S.rows("SELECT origin, COALESCE(verdict,'Untested') verdict, certainty, COUNT(*) n FROM claim "
                      "GROUP BY 1,2,3 ORDER BY 1,2,3")
    statements = S.rows("SELECT statement_id, short, verdict FROM consensus_statement")
    datasets = S.rows("SELECT dataset, cohort, modality, regions, n_donors, donors_by_disease FROM dataset")
    return dict(build=info[0] if info else {}, graph=g.stats(), claims_by_verdict=verdicts,
                consensus_statements=statements, datasets=datasets,
                cell_types=sorted(n.split(":", 1)[1] for n in g.nodes if n.startswith("cell:")),
                regions=sorted(n.split(":", 1)[1] for n in g.nodes if n.startswith("region:")),
                themes=sorted(n.split(":", 1)[1] for n in g.nodes if n.startswith("theme:")),
                relations={r: m for r, (_, _, m) in kg.RELATIONS.items()})


def _score(text, terms):
    t = text.lower()
    return sum(3 if re.search(rf"\b{re.escape(w)}\b", t) else 1 for w in terms if w in t)


def t_search(query, types=None, limit=20):
    limit = _int(limit, 20, 1, MAX_ROWS)
    terms = [w for w in re.findall(r"[\w\-]+", query.lower()) if len(w) > 1]
    if not terms:
        raise ToolError("Empty query.")
    g, hits = S.g, []
    for nid, n in g.nodes.items():
        if types and n["type"] not in types:
            continue
        text = " ".join([nid, str(n["label"])] + [str(n.get(k, "")) for k in
                                                   ("short", "subject", "cell_state", "title", "answer")])
        s = _score(text, terms)
        if s:
            hits.append((s + min(g.degree(nid), 50) / 100, nid))
    hits.sort(reverse=True)
    return dict(n_hits=len(hits), results=[_node_brief(n) for _, n in hits[:limit]])


def t_search_claims(query="", cell_type=None, region=None, gene=None, paper=None, verdict=None, origin=None,
                    theme=None, limit=25):
    limit = _int(limit, 25, 1, MAX_ROWS)
    where, args = [], []
    for col, val in (("cell_types", cell_type), ("regions", region), ("genes", gene), ("themes", theme)):
        if val:
            where.append(f"(';' || REPLACE({col}, '; ', ';') || ';') LIKE ?")
            args.append(f"%;{val.strip()};%")
    if paper:
        where.append("paper = ?"); args.append(paper)
    if verdict:
        where.append("COALESCE(verdict,'Untested') LIKE ?"); args.append(f"{verdict}%")
    if origin:
        where.append("origin = ?"); args.append(origin)
    rows = S.rows(f"SELECT {CLAIM_BRIEF}, genes, themes, context_text FROM claim"
                  + (" WHERE " + " AND ".join(where) if where else ""), args)
    terms = [w for w in re.findall(r"[\w\-]+", (query or "").lower()) if len(w) > 1]
    if terms:
        rows = [(r, _score(" ".join(str(r[k] or "") for k in ("claim_id", "statement", "genes", "themes",
                                                              "context_text", "cell_types", "regions")), terms))
                for r in rows]
        rows = [r for r, s in sorted(rows, key=lambda x: -x[1]) if s]
    return dict(n_matches=len(rows), claims=[_clean(r, drop=("context_text", "themes")) for r in rows[:limit]])


def t_get_claim(claim_id):
    claim_id = claim_id.removeprefix("claim:")
    r = S.rows("SELECT * FROM claim WHERE claim_id = ?", (claim_id,))
    if not r:
        raise ToolError(f"No claim {claim_id!r}; use search_claims.")
    c = _clean(r[0])
    tests = S.rows("SELECT dataset, role, outcome, effect, z, p, n_case, n_ref, region_match, regions_tested, note "
                   "FROM claim_dataset WHERE claim_id = ? AND outcome != 'not_testable' ORDER BY role, dataset",
                   (claim_id,))
    untestable = S.rows("SELECT dataset, note FROM claim_dataset WHERE claim_id = ? AND outcome = 'not_testable'",
                        (claim_id,))
    rel = S.rows("SELECT relation, claim_b target, basis, note FROM claim_relation WHERE claim_a = ? UNION ALL "
                 "SELECT relation || ' (from other claim)', claim_a, basis, note FROM claim_relation WHERE claim_b = ?",
                 (claim_id, claim_id))
    nid = f"claim:{claim_id}"
    bg = [dict(record=u, relation=r, statement=S.g.nodes[u]["label"]) for u, r, *_ in S.g.incident(nid)
          if S.g.nodes[u]["type"] == "background"]
    disc = S.rows("SELECT finding, evidence, found_by FROM source_discrepancy WHERE claim_id = ?", (claim_id,))
    return dict(claim=c, tests=[_clean(t) for t in tests],
                not_testable_in=[f"{u['dataset']}: {u['note']}" if u["note"] else u["dataset"] for u in untestable],
                relations=rel, reference_knowledge=bg, source_discrepancies=disc)


def t_get_dataset(dataset):
    dataset = dataset.removeprefix("dataset:")
    r = S.rows("SELECT * FROM dataset WHERE dataset = ? OR cohort = ?", (dataset, dataset))
    if not r:
        raise ToolError(f"No dataset {dataset!r}. Datasets: " + ", ".join(x["dataset"] for x in
                                                                         S.rows("SELECT dataset FROM dataset")))
    d = _clean(r[0])
    groups = S.rows("SELECT region, disease, n_donors FROM donor_group WHERE dataset = ?", (d["dataset"],))
    made = S.rows(f"SELECT {CLAIM_BRIEF} FROM claim WHERE source_dataset = ?", (d["dataset"],))
    tally = S.rows("SELECT role, outcome, COUNT(*) n FROM claim_dataset WHERE dataset = ? GROUP BY 1,2",
                   (d["dataset"],))
    return dict(dataset=d, donor_groups=groups, outcomes_of_tests_run_here=tally,
                claims_from_this_dataset=[_clean(x) for x in made])


def t_get_statement(statement_id):
    sid = statement_id.removeprefix("statement:")
    r = S.rows("SELECT * FROM consensus_statement WHERE statement_id = ?", (sid,))
    if not r:
        hits = S.g.resolve(statement_id, {"statement"})
        if not hits:
            raise ToolError("No such statement; see overview for the list.")
        sid = hits[0].split(":", 1)[1]
        r = S.rows("SELECT * FROM consensus_statement WHERE statement_id = ?", (sid,))
    links = S.rows(f"SELECT r.relation, r.note, {', '.join('c.' + x.strip() for x in CLAIM_BRIEF.split(','))} "
                   "FROM claim_relation r JOIN claim c ON c.claim_id = r.claim_a WHERE r.claim_b = ? "
                   "ORDER BY r.relation DESC, c.certainty", (sid,))
    bg = [dict(record=u, relation=rel, statement=S.g.nodes[u]["label"])
          for u, rel, *_ in S.g.incident(f"statement:{sid}") if S.g.nodes[u]["type"] == "background"]
    return dict(statement=_clean(r[0]), claims=[_clean(x) for x in links], reference_knowledge=bg)


def t_gene(symbol, include_signatures=True, limit=40):
    sym, limit = symbol.strip().upper(), _int(limit, 40, 1, MAX_ROWS)
    named = S.rows(f"SELECT {CLAIM_BRIEF} FROM claim WHERE (';' || REPLACE(genes,'; ',';') || ';') LIKE ?",
                   (f"%;{sym};%",))
    out = dict(gene=sym, claims_naming_gene=[_clean(x) for x in named])
    if include_signatures:
        sig = S.rows("SELECT cs.claim_id, s.signature_id, s.cell_type, s.cell_state, s.comparison, sg.rank, sg.stat, "
                     "sg.padj FROM signature_gene sg JOIN signature s ON s.signature_id = sg.signature_id "
                     "JOIN claim_signature cs ON cs.signature_id = sg.signature_id WHERE sg.gene = ? "
                     "ORDER BY sg.rank LIMIT ?", (sym, limit))
        out["claims_whose_signature_contains_gene"] = [_clean(x) for x in sig]
    if not named and not out.get("claims_whose_signature_contains_gene"):
        out["note"] = "Gene not named by any claim nor in any claim signature."
    return out


def t_neighbours(node, hops=1, relations=None, types=None, limit=80):
    nid = _resolve_one(node)
    hops, limit = _int(hops, 1, 1, 3, "hops"), _int(limit, 80, 1, 500)
    dist, edges = S.g.neighbourhood(nid, hops, set(relations) if relations else None,
                                    set(types) if types else None, limit=limit)
    return dict(center=_node_brief(nid),
                nodes=[dict(_node_brief(n), hop=h) for n, h in sorted(dist.items(), key=lambda x: (x[1], x[0]))
                       if n != nid],
                edges=[dict(source=s, relation=r, target=d, **{k: a[k] for k in ("outcome", "role", "note") if k in a})
                       for s, r, d, a in edges],
                note=f"{len(dist) - 1} neighbours within {hops} hop(s)" + (" (capped)" if len(dist) >= limit else ""))


def t_connect(a, b, max_hops=5, k=3):
    na, nb = _resolve_one(a), _resolve_one(b)
    max_hops, k = _int(max_hops, 5, 1, 6, "max_hops"), _int(k, 3, 1, 5, "k")
    ps = S.g.paths(na, nb, max_hops=max_hops, k=k)
    if not ps:
        return dict(a=na, b=nb, paths=[], note=f"No path within {max_hops} hops (themes are not used as links).")
    out = []
    for p in ps:
        cur, steps = na, []
        for s, r, d, attrs in p:
            nxt = d if s == cur else s
            steps.append(dict(**{"from": cur, "relation": r, "direction": "->" if s == cur else "<-", "to": nxt},
                              **{k: attrs[k] for k in ("outcome", "role") if k in attrs}))
            cur = nxt
        out.append(dict(steps=steps, nodes={n: S.g.nodes[n]["label"] for n in
                                            [na] + [st["to"] for st in steps]}))
    return dict(a=na, b=nb, paths=out)


def t_sql(query, limit=200):
    q, limit = query.strip().rstrip(";"), _int(limit, 200, 1, MAX_ROWS)
    if not re.match(r"(?is)^\s*(select|with|pragma\s+table_info)\b", q) or ";" in q:
        raise ToolError("Only a single SELECT / WITH statement (or PRAGMA table_info) is allowed.")
    try:
        with S.con(SQL_SECONDS) as c:
            cur = c.execute(q)
            cols = [d[0] for d in cur.description or []]
            rows = cur.fetchmany(limit + 1)
    except sqlite3.OperationalError as e:
        raise ToolError(f"SQL error: {e}" + (f" (queries are stopped after {SQL_SECONDS:g} s)" if "interrupt" in str(e)
                                            else "")) from None
    return dict(columns=cols, rows=[list(r) for r in rows[:limit]], truncated=len(rows) > limit)


def t_schema():
    with S.con() as c:
        tabs = [r[0] for r in c.execute("SELECT name FROM sqlite_master WHERE type='table' ORDER BY name")]
        return {t: [f"{r[1]} {r[2]}" for r in c.execute(f'PRAGMA table_info("{t}")')] for t in tabs}


def t_list_questions(status=None):
    rows = S.rows("SELECT id, text, asked_by, date, status, answer FROM question"
                  + (" WHERE status = ?" if status else ""), (status,) if status else ())
    pending = []
    if QUESTIONS.exists() and not READONLY:  # questions added since the last store build (local use only)
        import yaml
        ids = {r["id"] for r in rows}
        pending = [dict(id=x["id"], text=x["text"], status=x.get("status"), date=x.get("date"), note="not yet in store")
                   for x in (yaml.safe_load(QUESTIONS.read_text()) or []) if x.get("id") not in ids
                   and (not status or x.get("status") == status)]
    return dict(questions=rows + pending)


def t_add_question(text, asked_by="mcp client", context=""):
    if READONLY:
        raise ToolError(f"This server is read-only ({cfg.ENV_PREFIX}_MCP_READONLY=1).")
    import yaml
    text = " ".join(str(text).split())
    if len(text) < 8:
        raise ToolError("Question too short.")
    items = (yaml.safe_load(QUESTIONS.read_text()) or []) if QUESTIONS.exists() else []
    if any(x.get("text", "").lower() == text.lower() for x in items):
        raise ToolError("This question is already in the log.")
    slug = "-".join(re.findall(r"[a-z0-9]+", text.lower())[:6])
    qid = f"q-mcp-{date.today():%Y%m%d}-{slug}"
    entry = dict(id=qid, text=text, asked_by=f"mcp: {asked_by}"[:80], date=str(date.today()), status="open",
                 answer="", records=[], proposed={}, **({"context": context} if context else {}))
    items = sorted([x for x in items if x.get("id") != qid] + [entry], key=lambda x: x["id"])
    QUESTIONS.write_text(yaml.safe_dump(items, sort_keys=False, allow_unicode=True, width=120))
    return dict(added=entry, note="Logged in claims/questions.yaml; the atlas question loop (atlas-questions skill) "
                                  "picks it up for analysis.")


S_STR, S_INT = {"type": "string"}, {"type": "integer"}
S_LIST = {"type": "array", "items": {"type": "string"}}
TOOLS = {
    "overview": (t_overview, "What the store holds: build, counts, claims by verdict and certainty, consensus "
                 "statements, datasets, cell types, regions, themes and graph relation types. Call first.", {}, []),
    "search": (t_search, "Keyword search across every record in the knowledge graph (claims, findings, genes, cell "
               "types, regions, datasets, papers, consensus statements, reference knowledge, questions). Returns node "
               "ids for the other tools.",
               {"query": S_STR, "types": {**S_LIST, "description": f"restrict to node types: {kg.NODE_TYPES}"},
                "limit": S_INT}, ["query"]),
    "search_claims": (t_search_claims, "Find claims by keyword and/or filters; returns verdict, certainty and pooled "
                      "effect (est, 95% CI lo-hi, approx. SMD, positive = as claimed).",
                      {"query": S_STR, "cell_type": {**S_STR, "description": f"e.g. {EX['cell_types']}"},
                       "region": {**S_STR, "description": f"e.g. {EX['regions']}"}, "gene": S_STR,
                       "paper": {**S_STR, "description": f"paper id, e.g. {EX['paper']}; this-study for findings"},
                       "verdict": {**S_STR, "description": "Supported, Partly supported, Inconclusive, No effect, "
                                   "Contradicted, Untested"},
                       "origin": {**S_STR, "description": "published or finding"}, "theme": S_STR, "limit": S_INT},
                      []),
    "get_claim": (t_get_claim, "Everything about one claim: statement, quote and location in the paper, full context "
                  "(organism, region, cell type, contrast, cohort), rating with reasons, per-dataset re-test results, "
                  "relations, reference knowledge and known source discrepancies.",
                  {"claim_id": {**S_STR, "description": f"e.g. {EX['claim']} or {EX['finding']}"}}, ["claim_id"]),
    "get_statement": (t_get_statement, "A consensus statement with its verdict, the claims that support or oppose it "
                      "(with their ratings) and the reference knowledge behind it.",
                      {"statement_id": {**S_STR, "description": f"e.g. {EX['statement']}, or words of the statement"}},
                      ["statement_id"]),
    "get_dataset": (t_get_dataset, "One dataset/cohort: assay, regions, donors by disease, ancestry, brain banks, "
                    "staging, sampling, the claims made from it and the outcomes of tests run in it.",
                    {"dataset": S_STR}, ["dataset"]),
    "gene": (t_gene, "Claims that name a gene, and claims whose published signature contains it (with rank).",
             {"symbol": S_STR, "include_signatures": {"type": "boolean"}, "limit": S_INT}, ["symbol"]),
    "neighbours": (t_neighbours, "Knowledge-graph neighbourhood of a node (1-3 hops), optionally restricted to "
                   f"relation or node types. Node by id (claim:{EX['claim']}, gene:{EX['gene']}, cell:{EX['cell']}, "
                   f"region:{EX['region']}, statement:{EX['statement']}, dataset:{EX['dataset']}, "
                   f"background:{EX['background']}) or by name.",
                   {"node": S_STR, "hops": S_INT,
                    "relations": {**S_LIST, "description": f"relation types: {list(kg.RELATIONS)}"},
                    "types": {**S_LIST, "description": f"node types: {kg.NODE_TYPES}"}, "limit": S_INT}, ["node"]),
    "connect": (t_connect, "How two records are linked: the most specific paths between two nodes in the knowledge "
                "graph (hub nodes cost more, themes are skipped).",
                {"a": S_STR, "b": S_STR, "max_hops": S_INT, "k": S_INT}, ["a", "b"]),
    "schema": (t_schema, "Tables and columns of the SQLite claim store, for use with `sql`.", {}, []),
    "sql": (t_sql, "Read-only SQL (SQLite) over the claim store. Main tables: claim, claim_dataset (re-tests), "
            "claim_relation, consensus_statement, dataset, donor_group, source, signature, signature_gene, "
            "claim_signature, background, question, vocab_*. Call `schema` for columns.",
            {"query": S_STR, "limit": S_INT}, ["query"]),
    "list_questions": (t_list_questions, "Questions asked of the atlas and their status and answers.",
                       {"status": {**S_STR, "description": "open, answered, analysed, needs_analysis, no_data"}}, []),
    "add_question": (t_add_question, "Log a question the store cannot answer yet, for a new donor-level analysis "
                     "(writes claims/questions.yaml). Use only when the user wants it logged.",
                     {"text": S_STR, "asked_by": S_STR, "context": S_STR}, ["text"]),
}
READ_ONLY_TOOLS = set(TOOLS) - {"add_question"}


def tool_list():
    out = []
    for name, (_, desc, props, req) in TOOLS.items():
        if READONLY and name not in READ_ONLY_TOOLS:
            continue
        out.append(dict(name=name, description=desc,
                        inputSchema={"type": "object", "properties": props, "required": req},
                        annotations={"readOnlyHint": name in READ_ONLY_TOOLS, "openWorldHint": False}))
    return out


def call_tool(name, args):
    if name not in TOOLS or (READONLY and name not in READ_ONLY_TOOLS):
        raise ToolError(f"Unknown tool {name!r}.")
    fn = TOOLS[name][0]
    if args is not None and not isinstance(args, dict):
        raise ToolError("arguments must be a JSON object")
    try:
        res = fn(**(args or {}))
    except TypeError as e:
        raise ToolError(f"Bad arguments for {name}: {e}") from e
    return json.dumps(res, default=str, ensure_ascii=False, separators=(",", ":"))


# ------------------------------------------------------------------ JSON-RPC over stdio
def handle(msg):
    mid, method, params = msg.get("id"), msg.get("method"), msg.get("params") or {}
    if method is None or mid is None:  # response or notification (initialized, cancelled): nothing to send
        return None
    try:
        if method == "initialize":
            v = params.get("protocolVersion")
            result = dict(protocolVersion=v if v in PROTOCOLS else PROTOCOLS[0],
                          capabilities={"tools": {"listChanged": False}},
                          serverInfo={"name": ABOUT["server_name"], "title": ABOUT["name"], "version": VERSION},
                          instructions=INSTRUCTIONS)
        elif method == "ping":
            result = {}
        elif method == "tools/list":
            result = {"tools": tool_list()}
        elif method == "tools/call":
            try:
                result = {"content": [{"type": "text", "text": call_tool(params.get("name"), params.get("arguments"))}],
                          "isError": False}
            except ToolError as e:
                result = {"content": [{"type": "text", "text": str(e)}], "isError": True}
            except Exception as e:  # noqa: BLE001 - report to the client instead of dying
                traceback.print_exc(file=sys.stderr)
                result = {"content": [{"type": "text", "text": f"{type(e).__name__}: {e}"}], "isError": True}
        elif method in ("resources/list", "prompts/list", "resources/templates/list"):
            result = {method.split("/")[0]: [], **({"resourceTemplates": []} if "templates" in method else {})}
        else:
            return {"jsonrpc": "2.0", "id": mid, "error": {"code": -32601, "message": f"Method not found: {method}"}}
        return {"jsonrpc": "2.0", "id": mid, "result": result}
    except Exception as e:  # noqa: BLE001
        traceback.print_exc(file=sys.stderr)
        return {"jsonrpc": "2.0", "id": mid, "error": {"code": -32603, "message": str(e)}}


def main():
    print(f"{ABOUT['server_name']} MCP server {VERSION}: {DB}", file=sys.stderr)
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        try:
            msg = json.loads(line)
        except ValueError:
            out = {"jsonrpc": "2.0", "id": None, "error": {"code": -32700, "message": "Parse error"}}
        else:
            out = [r for r in map(handle, msg) if r] if isinstance(msg, list) else handle(msg)
        if out:
            sys.stdout.write(json.dumps(out, default=str, ensure_ascii=False) + "\n")
            sys.stdout.flush()


if __name__ == "__main__":
    main()

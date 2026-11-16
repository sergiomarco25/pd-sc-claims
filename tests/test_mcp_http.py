"""End-to-end test of the remote server (src/pdrep/mcp_http.py): MCP over HTTP, REST, limits and safety.

Starts the server on a free local port against a claim store, runs every check, stops it. Standard library only.

    python3 -I tests/test_mcp_http.py [--db path/to/claim_store.sqlite] [--url http://host:port]
    pytest tests/test_mcp_http.py          (uses $PDREP_STORE_DB, else the test fixture, else the configured store)

With --url it tests an already running server (for example the deployed one) instead of starting its own.
"""
import argparse
import json
import os
import socket
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "tests" / "fixture" / "claim_store.sqlite"
PROTO = "2025-06-18"
N_TOOLS = 12  # read-only tool set (add_question hidden)
CLAIM = "xu2023-C11"  # present in the full store and in the fixture


def default_db() -> Path:
    if os.environ.get("PDREP_STORE_DB"):
        return Path(os.environ["PDREP_STORE_DB"])
    if FIXTURE.is_file():
        return FIXTURE
    sys.path.insert(0, str(ROOT / "src"))
    from pdrep import config as cfg
    return cfg.STORE_DB


class Client:
    def __init__(self, base, token=""):
        self.base, self.token = base.rstrip("/"), token

    def __call__(self, method, path, body=None, headers=None):
        h = {"Content-Type": "application/json", "Accept": "application/json, text/event-stream", **(headers or {})}
        if self.token:
            h["Authorization"] = f"Bearer {self.token}"
        data = None if body is None else json.dumps(body).encode()
        r = urllib.request.Request(self.base + path, data=data, method=method, headers=h)
        try:
            with urllib.request.urlopen(r, timeout=60) as f:
                raw = f.read()
                is_json = raw and "json" in f.headers.get("Content-Type", "")
                return f.status, (json.loads(raw) if is_json else raw.decode()), dict(f.headers)
        except urllib.error.HTTPError as e:
            raw = e.read()
            return e.code, (json.loads(raw) if raw else None), dict(e.headers)


def text(r):
    return json.loads(r["result"]["content"][0]["text"])


def checks(req):
    """Yield (label, passed, detail)."""
    s, r, _ = req("GET", "/health")
    yield "health", s == 200 and r["status"] == "ok", ""
    s, r, _ = req("GET", "/")
    yield "landing text", s == 200 and "/mcp" in r, ""
    init = {"jsonrpc": "2.0", "id": 1, "method": "initialize",
            "params": {"protocolVersion": PROTO, "capabilities": {}, "clientInfo": {"name": "test", "version": "0"}}}
    s, r, _ = req("POST", "/mcp", init)
    yield "mcp initialize", s == 200 and r["result"]["protocolVersion"] == PROTO and "instructions" in r["result"], ""
    hdr = {"MCP-Protocol-Version": PROTO}
    s, r, _ = req("POST", "/mcp", {"jsonrpc": "2.0", "method": "notifications/initialized"}, hdr)
    yield "notification -> 202", s == 202, s
    s, r, _ = req("POST", "/mcp", {"jsonrpc": "2.0", "id": 2, "method": "tools/list"}, hdr)
    names = [t["name"] for t in r["result"]["tools"]]
    yield "tools/list is read-only", s == 200 and "add_question" not in names and len(names) == N_TOOLS, names
    call = {"jsonrpc": "2.0", "id": 3, "method": "tools/call", "params": {"name": "get_claim", "arguments": {"claim_id": CLAIM}}}
    s, r, _ = req("POST", "/mcp", call)
    yield "tools/call get_claim", s == 200 and not r["result"]["isError"] and text(r)["claim"]["claim_id"] == CLAIM, ""
    call = {"jsonrpc": "2.0", "id": 4, "method": "tools/call", "params": {"name": "add_question", "arguments": {"text": "x"}}}
    s, r, _ = req("POST", "/mcp", call)
    yield "add_question refused", r["result"]["isError"], ""
    batch = [{"jsonrpc": "2.0", "id": 5, "method": "ping"},
             {"jsonrpc": "2.0", "id": 6, "method": "tools/call",
              "params": {"name": "search", "arguments": {"query": "microglia", "limit": 5}}}]
    s, r, _ = req("POST", "/mcp", batch)
    yield "batch", s == 200 and isinstance(r, list) and len(r) == 2 and not r[1]["result"]["isError"], ""
    s, r, _ = req("POST", "/mcp", init, {"MCP-Protocol-Version": "1999-01-01"})
    yield "bad protocol header -> 400", s == 400, s
    s, r, _ = req("GET", "/mcp")
    yield "GET /mcp -> 405", s == 405, s
    q = "WITH RECURSIVE x(n) AS (SELECT 1 UNION ALL SELECT n+1 FROM x) SELECT count(*) FROM x"
    s, r, _ = req("POST", "/mcp", {"jsonrpc": "2.0", "id": 7, "method": "tools/call",
                                   "params": {"name": "sql", "arguments": {"query": q}}})
    yield "runaway sql stopped", r["result"]["isError"] and "stopped" in r["result"]["content"][0]["text"], ""
    s, r, _ = req("POST", "/mcp", {"jsonrpc": "2.0", "id": 8, "method": "tools/call",
                                   "params": {"name": "sql", "arguments": {"query": "DELETE FROM claim"}}})
    yield "write sql refused", r["result"]["isError"], ""
    s, r, _ = req("POST", "/mcp", {"jsonrpc": "2.0", "id": 9, "method": "tools/call",
                                   "params": {"name": "search", "arguments": {"query": "microglia", "limit": "100000"}}})
    yield "limit capped", len(text(r)["results"]) <= 1000, ""
    s, r, _ = req("GET", "/api/tools")
    yield "REST tool list", s == 200 and len(r["tools"]) == N_TOOLS, ""
    s, r, _ = req("POST", "/api/tools/search_claims", {"cell_type": "Microglia"})
    yield "REST search_claims", s == 200 and "claims" in r, ""
    s, r, _ = req("POST", "/api/tools/get_claim", {"claim_id": "nope-C99"})
    yield "REST error -> 400", s == 400 and "error" in r, s
    s, r, _ = req("POST", "/api/tools/nosuch", {})
    yield "REST unknown tool -> 404", s == 404, s
    s, r, _ = req("POST", "/api/tools/gene", ["not", "object"])
    yield "REST non-object body -> 400", s == 400, s
    s, r, _ = req("GET", "/openapi.json")
    yield "openapi", s == 200 and r["openapi"] == "3.1.0" and len(r["paths"]) == N_TOOLS, ""


def free_port():
    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        return s.getsockname()[1]


def start_server(db: Path):
    port = free_port()
    env = {**os.environ, "PDREP_STORE_DB": str(db), "PDREP_HTTP_HOST": "127.0.0.1", "PORT": str(port),
           "PDREP_RATE_PER_MIN": "100000"}
    env.pop("PDREP_MCP_ALLOW_WRITE", None)
    env.pop("PDREP_MCP_TOKEN", None)
    p = subprocess.Popen([sys.executable, "-I", str(ROOT / "src" / "pdrep" / "mcp_http.py")], env=env,
                         stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    base = f"http://127.0.0.1:{port}"
    for _ in range(100):
        try:
            urllib.request.urlopen(base + "/health", timeout=1)
            return p, base
        except OSError:
            if p.poll() is not None:
                raise RuntimeError("server exited: " + p.stderr.read().decode()[-2000:])
            time.sleep(0.2)
    p.kill()
    raise RuntimeError("server did not start")


def run(db=None, url=None, token=""):
    proc = None
    if not url:
        proc, url = start_server(Path(db or default_db()))
    try:
        results = list(checks(Client(url, token)))
    finally:
        if proc:
            proc.terminate()
            proc.wait(timeout=10)
    return results


def test_http_server():  # pytest entry point
    failed = [(label, detail) for label, ok, detail in run() if not ok]
    assert not failed, failed


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--db")
    ap.add_argument("--url")
    ap.add_argument("--token", default="")
    a = ap.parse_args()
    res = run(a.db, a.url, a.token)
    for label, ok, detail in res:
        print(("PASS " if ok else "FAIL ") + label, detail if not ok else "")
    n_fail = sum(not ok for _, ok, _ in res)
    print(f"{len(res) - n_fail}/{len(res)} passed")
    sys.exit(1 if n_fail else 0)

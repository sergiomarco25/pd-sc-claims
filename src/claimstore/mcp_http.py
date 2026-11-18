"""Remote (HTTP) server for a claim store: the same tools as mcp_server.py, for any AI system.

Two interfaces on one port, standard library only (Python 3.9+, no install):

  POST /mcp                  Model Context Protocol, Streamable HTTP transport (JSON responses, stateless).
                             For MCP clients: Claude (claude.ai / Desktop / Code custom connectors), ChatGPT
                             connectors, the OpenAI Responses API `mcp` tool, Gemini CLI, Cursor, VS Code,
                             LangChain / LlamaIndex MCP adapters, ...
  GET  /openapi.json         OpenAPI 3.1 description of the REST interface below, for systems that call plain
                             HTTP functions (custom GPT Actions, function calling in any LLM API, scripts).
  GET  /api/tools            tool list with JSON Schemas
  POST /api/tools/<name>     call a tool; body = JSON object of arguments; response = the tool's JSON result
  GET  /health, GET /        liveness, short usage text

Read-only by default: add_question (the one tool that writes) is only offered when <PREFIX>_MCP_ALLOW_WRITE=1.
The store is opened read-only; sql queries are stopped after <PREFIX>_SQL_SECONDS (default 5 s).

Run locally:
    python3 -I src/claimstore/mcp_http.py                # http://127.0.0.1:8000/mcp
Run in a container / cloud service (see deploy/mcp/README.md):
    <PREFIX>_HTTP_HOST=0.0.0.0 PORT=8080 <PREFIX>_STORE_DB=/app/claim_store.sqlite python3 -I mcp_http.py

Environment (<PREFIX> = [project] env_prefix of the project config):
  <PREFIX>_STORE_DB          claim_store.sqlite to serve
  <PREFIX>_HTTP_HOST, PORT   bind address (default 127.0.0.1:8000; PORT is what Cloud Run / Render / Fly set)
  <PREFIX>_PUBLIC_URL        public base URL, written into openapi.json (e.g. https://claims.example.org)
  <PREFIX>_MCP_TOKEN         optional shared secret: if set, every /mcp and /api call needs "Authorization: Bearer <it>"
                          (scripts, OpenAI API, Cursor; claude.ai connectors use OAuth or no auth, so leave it unset
                          for a public read-only claude.ai connector)
  <PREFIX>_ALLOWED_ORIGINS   comma-separated browser origins allowed to call (default: any origin when bound publicly,
                          only localhost origins when bound to localhost, which blocks DNS-rebinding attacks)
  <PREFIX>_RATE_PER_MIN      requests per minute per client IP (default 120)
  <PREFIX>_MCP_ALLOW_WRITE=1 offer add_question (writes claims/questions.yaml next to the code; local use only)
"""
import hmac
import json
import os
import sys
import threading
import time
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlsplit

if __package__ in (None, ""):  # run as a file: make the package importable (repository src/, or /app in a container)
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from claimstore import config as cfg  # noqa: E402

if cfg.env("MCP_ALLOW_WRITE") != "1":
    os.environ[f"{cfg.ENV_PREFIX}_MCP_READONLY"] = "1"  # before import: mcp_server reads it at import time
from claimstore import mcp_server as M  # noqa: E402

HOST = cfg.env("HTTP_HOST", "127.0.0.1")
PORT = int(os.environ.get("PORT", cfg.env("HTTP_PORT", "8000")))
PUBLIC_URL = (cfg.env("PUBLIC_URL", "") or "").rstrip("/")
TOKEN = cfg.env("MCP_TOKEN", "") or ""
ORIGINS = {o.strip().rstrip("/") for o in (cfg.env("ALLOWED_ORIGINS", "") or "").split(",") if o.strip()}
RATE = int(cfg.env("RATE_PER_MIN", "120"))
MAX_BODY = 1 << 20  # 1 MiB
LOCAL = HOST in ("127.0.0.1", "localhost", "::1")


# ------------------------------------------------------------------ OpenAPI
def openapi():
    paths = {}
    for t in M.tool_list():
        paths[f"/api/tools/{t['name']}"] = {"post": {
            "operationId": t["name"], "summary": t["description"][:300], "description": t["description"],
            "requestBody": {"required": bool(t["inputSchema"]["required"]),
                            "content": {"application/json": {"schema": t["inputSchema"]}}},
            "responses": {"200": {"description": "Tool result (JSON)",
                                  "content": {"application/json": {"schema": {"type": "object"}}}},
                          "400": {"description": "Bad arguments or no match",
                                  "content": {"application/json": {"schema": {
                                      "type": "object", "properties": {"error": {"type": "string"}}}}}}}}}
    doc = {"openapi": "3.1.0",
           "info": {"title": M.ABOUT["name"], "version": M.VERSION, "description": M.INSTRUCTIONS},
           "paths": paths}
    if PUBLIC_URL:
        doc["servers"] = [{"url": PUBLIC_URL}]
    if TOKEN:
        doc["components"] = {"securitySchemes": {"bearer": {"type": "http", "scheme": "bearer"}}}
        doc["security"] = [{"bearer": []}]
    return doc


HELP = """{name} - MCP and REST server {v}

MCP (Streamable HTTP):  POST {b}/mcp
REST tools:             GET  {b}/api/tools ; POST {b}/api/tools/<name> with a JSON body of arguments
OpenAPI description:    GET  {b}/openapi.json
Tools: {tools}
"""


# ------------------------------------------------------------------ rate limit (per client IP, token bucket)
_buckets, _blk = {}, threading.Lock()


def allow(ip):
    now = time.monotonic()
    with _blk:
        tokens, last = _buckets.get(ip, (RATE, now))
        tokens = min(RATE, tokens + (now - last) * RATE / 60)
        if tokens < 1:
            _buckets[ip] = (tokens, now)
            return False
        _buckets[ip] = (tokens - 1, now)
        if len(_buckets) > 10000:  # forget idle clients
            for k in [k for k, (_, t) in _buckets.items() if now - t > 600]:
                _buckets.pop(k, None)
        return True


# ------------------------------------------------------------------ HTTP
class Handler(BaseHTTPRequestHandler):
    server_version = f"{M.ABOUT['server_name']}/{M.VERSION}"
    protocol_version = "HTTP/1.1"

    def log_message(self, fmt, *args):  # one line per request on stderr, no bodies
        sys.stderr.write(f"{self.client_ip()} {self.command} {self.path.split('?')[0]} {fmt % args}\n")

    def client_ip(self):
        fwd = self.headers.get("X-Forwarded-For", "")  # behind a cloud load balancer
        return fwd.split(",")[0].strip() if fwd else self.client_address[0]

    # -- responses
    def cors(self):
        o = self.headers.get("Origin")
        if o and self.origin_ok():
            self.send_header("Access-Control-Allow-Origin", o)
            self.send_header("Vary", "Origin")
            self.send_header("Access-Control-Allow-Headers",
                             "Content-Type, Authorization, Accept, Mcp-Session-Id, MCP-Protocol-Version")
            self.send_header("Access-Control-Allow-Methods", "GET, POST, OPTIONS")
            self.send_header("Access-Control-Expose-Headers", "Mcp-Session-Id")

    def send(self, code, body=None, ctype="application/json", extra=()):
        data = b"" if body is None else (body if isinstance(body, bytes) else
                                         (json.dumps(body, ensure_ascii=False, default=str) if ctype ==
                                          "application/json" else body).encode())
        self.send_response(code)
        if body is not None:
            self.send_header("Content-Type", ctype + ("; charset=utf-8" if not ctype.endswith("json") else ""))
        self.send_header("Content-Length", str(len(data)))
        self.send_header("Cache-Control", "no-store")
        for k, v in extra:
            self.send_header(k, v)
        self.cors()
        self.end_headers()
        if data and self.command != "HEAD":
            self.wfile.write(data)

    # -- checks
    def origin_ok(self):
        o = self.headers.get("Origin")
        if not o:
            return True  # not a browser
        o = o.rstrip("/")
        if ORIGINS:
            return o in ORIGINS or "*" in ORIGINS
        if LOCAL:  # DNS-rebinding protection for a server on localhost
            return urlsplit(o).hostname in ("localhost", "127.0.0.1", "::1")
        return True

    def guard(self):
        """Origin, auth and rate checks shared by /mcp and /api; returns False after sending an error."""
        if not self.origin_ok():
            self.send(403, {"error": "origin not allowed"})
            return False
        if TOKEN:
            got = self.headers.get("Authorization", "")
            if not hmac.compare_digest(got.encode(), f"Bearer {TOKEN}".encode()):
                self.send(401, {"error": "missing or wrong bearer token"}, extra=[("WWW-Authenticate", "Bearer")])
                return False
        if not allow(self.client_ip()):
            self.send(429, {"error": "too many requests"}, extra=[("Retry-After", "10")])
            return False
        return True

    def body(self):
        n = int(self.headers.get("Content-Length") or 0)
        if n > MAX_BODY:
            raise ValueError("request too large")
        raw = self.rfile.read(n) if n else b""
        return json.loads(raw) if raw.strip() else None

    def base(self):
        if PUBLIC_URL:
            return PUBLIC_URL
        proto = self.headers.get("X-Forwarded-Proto", "http")
        return f"{proto}://{self.headers.get('Host', f'{HOST}:{PORT}')}"

    # -- verbs
    def do_OPTIONS(self):
        self.send(204 if self.origin_ok() else 403)

    def do_HEAD(self):
        self.do_GET()

    def do_GET(self):
        path = urlsplit(self.path).path.rstrip("/") or "/"
        if path == "/health":
            return self.send(200, {"status": "ok", "version": M.VERSION})
        if path == "/":
            return self.send(200, HELP.format(name=M.ABOUT['name'], v=M.VERSION, b=self.base(),
                                              tools=", ".join(t["name"] for t in M.tool_list())), "text/plain")
        if path == "/openapi.json":
            return self.send(200, openapi())
        if path == "/mcp":  # no server-initiated stream: the spec allows 405 here
            return self.send(405, {"error": "use POST for MCP messages"}, extra=[("Allow", "POST, OPTIONS")])
        if path == "/api/tools":
            if not self.guard():
                return
            return self.send(200, {"tools": M.tool_list(), "instructions": M.INSTRUCTIONS})
        self.send(404, {"error": "not found"})

    def do_DELETE(self):
        self.send(405, {"error": "stateless server: no sessions to delete"}, extra=[("Allow", "POST, OPTIONS")])

    def do_POST(self):
        path = urlsplit(self.path).path.rstrip("/")
        if path not in ("/mcp",) and not path.startswith("/api/tools/"):
            return self.send(404, {"error": "not found"})
        if not self.guard():
            return
        try:
            msg = self.body()
        except ValueError as e:
            if path == "/mcp":
                return self.send(400, {"jsonrpc": "2.0", "id": None, "error": {"code": -32700, "message": str(e)}})
            return self.send(400, {"error": f"invalid JSON: {e}"})
        if path == "/mcp":
            return self.mcp(msg)
        name = path.rsplit("/", 1)[1]
        try:
            return self.send(200, json.loads(M.call_tool(name, msg or {})))
        except M.ToolError as e:
            return self.send(404 if str(e).startswith("Unknown tool") else 400, {"error": str(e)})
        except Exception as e:  # noqa: BLE001
            M.traceback.print_exc(file=sys.stderr)
            return self.send(500, {"error": f"{type(e).__name__}: {e}"})

    def mcp(self, msg):
        v = self.headers.get("MCP-Protocol-Version")
        if v and v not in M.PROTOCOLS:
            return self.send(400, {"jsonrpc": "2.0", "id": None,
                                   "error": {"code": -32600, "message": f"unsupported protocol version {v}"}})
        if not isinstance(msg, (dict, list)) or (isinstance(msg, list) and not msg):
            return self.send(400, {"jsonrpc": "2.0", "id": None, "error": {"code": -32600, "message": "invalid request"}})
        out = [r for r in map(M.handle, msg) if r] if isinstance(msg, list) else M.handle(msg)
        if not out:  # only notifications / responses
            return self.send(202)
        self.send(200, out)


def main():
    srv = ThreadingHTTPServer((HOST, PORT), Handler)
    srv.daemon_threads = True
    M.S.g  # build the graph once before serving
    print(f"{M.ABOUT['server_name']} {M.VERSION} on http://{HOST}:{PORT} (MCP at /mcp, REST at /api/tools, /openapi.json); "
          f"store {M.DB}; {'read-only' if M.READONLY else 'WRITES ENABLED'}; "
          f"{'token required' if TOKEN else 'no auth'}", file=sys.stderr)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    main()

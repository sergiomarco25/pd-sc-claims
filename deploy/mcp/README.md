# Remote MCP server for the PD claim store

One small web service gives any AI system read access to the claim store.

- **MCP** (Model Context Protocol, Streamable HTTP) at `/mcp`: the open standard used by Claude, ChatGPT, the OpenAI API, Gemini CLI, Cursor, VS Code and the LangChain / LlamaIndex adapters.
- **REST + OpenAPI** at `/api/tools/<tool>` and `/openapi.json`: for anything that calls plain HTTP functions, such as custom GPT Actions, function calling in any LLM API, R or Python scripts, and notebooks.

It offers the same 12 read-only tools as the local server (`src/pdrep/mcp_server.py`):

| Tool | What it does |
|---|---|
| `overview` | What the store holds |
| `search`, `search_claims` | Find records and claims |
| `get_claim`, `get_statement`, `get_dataset`, `gene` | Read one record |
| `neighbours`, `connect` | Explore the knowledge graph |
| `schema`, `sql` | Read-only SQL |
| `list_questions` | The question log |

`add_question` is switched off in the remote server.

**What it serves:** derived knowledge only. That means published statements and short quotes, context, ratings, per-dataset effect sizes and z / p values, and donor counts per group. There are no cells, expression matrices or donor-level records. The image holds a snapshot of `claim_store.sqlite` (about 19 MB) next to the code, which is standard-library Python with no dependencies.

## 1. Build the image

```bash
deploy/mcp/stage.sh                      # code + consistent snapshot of the store -> deploy/mcp/build/
docker build -t pd-claim-store deploy/mcp/build
docker run --rm -p 8080:8080 pd-claim-store   # try it: http://localhost:8080/
```

The cluster can't build images (rootless podman has no subuid ranges), so build on a laptop or in the cloud (below). Rerun `stage.sh` and redeploy after every store rebuild (`51b_claim_store.py`).

## 2. Put it online (pick one)

The service must be reachable over **HTTPS** from the internet, because claude.ai and ChatGPT connect from their own servers. The Helmholtz cluster is not reachable from outside.

**Google Cloud Run** (scales to zero, free tier is usually enough):

```bash
cd deploy/mcp/build
gcloud run deploy pd-claim-store --source . --region europe-west3 --allow-unauthenticated \
  --memory 512Mi --max-instances 3 \
  --set-env-vars PDREP_PUBLIC_URL=https://<the URL Cloud Run prints>
```

**Fly.io / Render / Railway:** deploy `deploy/mcp/build/` as a Dockerfile service with internal port 8080. They set `PORT` and provide HTTPS.

**A VM** (for example Helmholtz / de.NBI cloud):
1. Run `docker run -d --restart always -p 127.0.0.1:8080:8080 pd-claim-store`.
2. Put Caddy in front for automatic HTTPS. A single Caddyfile line is enough: `pd-claims.example.org { reverse_proxy 127.0.0.1:8080 }`.

Check it:

```bash
curl https://<host>/health
```

## 3. Connect AI systems

Replace `https://<host>` with the public URL. The MCP endpoint is `https://<host>/mcp`.

| System | How |
|---|---|
| **Claude** (claude.ai, Desktop, mobile) | Settings → Connectors → *Add custom connector* → URL `https://<host>/mcp`. On Team/Enterprise plans an Owner adds it once in the organisation settings, then members enable it. |
| **Claude Code** | `claude mcp add --transport http pd-claim-store https://<host>/mcp` |
| **ChatGPT** | Settings → Apps & Connectors → developer mode → create a connector with URL `https://<host>/mcp` (plan-dependent). |
| **OpenAI API** (Responses) | `tools=[{"type": "mcp", "server_label": "pd_claims", "server_url": "https://<host>/mcp", "require_approval": "never"}]` |
| **Gemini CLI** | `~/.gemini/settings.json`: `{"mcpServers": {"pd-claims": {"httpUrl": "https://<host>/mcp"}}}` |
| **Cursor** | `.cursor/mcp.json`: `{"mcpServers": {"pd-claims": {"url": "https://<host>/mcp"}}}` |
| **VS Code** (Copilot agent) | `.vscode/mcp.json`: `{"servers": {"pd-claims": {"type": "http", "url": "https://<host>/mcp"}}}` |
| **LangChain / LangGraph** | `langchain-mcp-adapters`: `MultiServerMCPClient({"pd": {"transport": "streamable_http", "url": "https://<host>/mcp"}})` |
| **Custom GPT Actions**, other function-calling LLMs | Import `https://<host>/openapi.json`. Each tool becomes one operation. |
| **Scripts** (Python, R, curl) | `curl -X POST https://<host>/api/tools/search_claims -H 'Content-Type: application/json' -d '{"cell_type":"Microglia","verdict":"Contradicted"}'` |

These client settings change often, so if one doesn't work, check that client's current documentation. Any client that supports "remote MCP over Streamable HTTP" will work.

## 4. Access and safety

- **Default: public and read-only, no sign-in.** This suits a resource meant for the field.
  - The store is opened read-only, and `sql` accepts only single SELECT statements, stopped after 5 s.
  - Results are capped at 1,000 rows.
  - Requests are limited to 120 per minute per IP.
  - Request bodies are limited to 1 MiB.
- **Shared secret:** set `PDREP_MCP_TOKEN=<secret>`, and every call then needs `Authorization: Bearer <secret>`. This works for scripts, the OpenAI API, Cursor, VS Code and Claude Code. The claude.ai and ChatGPT connector screens expect OAuth instead, so for them, keep the server public or add an OAuth proxy.
- **Browser origins:** `PDREP_ALLOWED_ORIGINS=https://a.org,https://b.org` restricts which web pages may call it from a browser.

| Setting | Default | Meaning |
|---|---|---|
| `PDREP_STORE_DB` | `/app/claim_store.sqlite` in the image | Store to serve |
| `PORT`, `PDREP_HTTP_HOST` | `8080`, `0.0.0.0` in the image | Bind address |
| `PDREP_PUBLIC_URL` | unset | Written into `openapi.json` |
| `PDREP_MCP_TOKEN` | unset | Optional bearer token |
| `PDREP_ALLOWED_ORIGINS` | any (localhost only when bound to localhost) | Browser origins |
| `PDREP_RATE_PER_MIN` | 120 | Per-IP rate limit |
| `PDREP_SQL_SECONDS` | 5 | Time limit per `sql` query |

## Run it locally without a container

```bash
python3 -I src/pdrep/mcp_http.py            # http://127.0.0.1:8000/mcp and /api/tools
```

## Tested

On 2026-11-13, against the staged `build/` contents run as in the image:
- MCP initialize, `tools/list` and `tools/call`, notifications, batches and the protocol-version header.
- REST calls, error codes and `openapi.json`.
- The SQL time limit, the row cap, origin checks, the token and the rate limit.
- The official MCP Python SDK client (`mcp` 2.3.0, Streamable HTTP), which negotiated protocol 2025-11-25 and called the tools.

Not yet tested: the container build itself (needs Docker off the cluster) and the hosted claude.ai and ChatGPT connectors.

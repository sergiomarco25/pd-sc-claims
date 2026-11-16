---
title: PD single-cell claim store
emoji: 🧠
colorFrom: blue
colorTo: indigo
sdk: docker
app_port: 8080
pinned: false
short_description: Rated PD single-cell claims for AI systems (MCP)
---

# PD single-cell claim store: MCP and REST server

This server holds the published conclusions of human post-mortem Parkinson's disease single-cell and single-nucleus studies, extracted as structured claims. Each claim was re-tested at donor level in its source dataset and in every other harmonised PD dataset, then rated with a verdict and a certainty. The server also holds consensus statements and reference knowledge.

It serves derived knowledge only: statements, context, ratings, per-dataset effect sizes and donor counts. There are no cells, expression matrices or donor-level records.

## Use it from an AI system

- **MCP** (Claude, ChatGPT, OpenAI API, Gemini CLI, Cursor, VS Code, LangChain and others): add `https://<this-space>.hf.space/mcp` as a remote MCP server or custom connector.
- **REST / OpenAPI** (custom GPT Actions, function calling, scripts): `https://<this-space>.hf.space/openapi.json`, `GET /api/tools`, `POST /api/tools/<tool>`.

The server offers 12 read-only tools:

| Tool | What it does |
|---|---|
| `overview` | What the store holds |
| `search`, `search_claims` | Find records and claims |
| `get_claim`, `get_statement`, `get_dataset`, `gene` | Read one record |
| `neighbours`, `connect` | Explore the knowledge graph |
| `schema`, `sql` | Read-only SQL |
| `list_questions` | The question log |

**Cold starts:** free Spaces sleep after about 2 days without traffic, so the first request after a pause can take about 30 s.

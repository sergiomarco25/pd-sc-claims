# Claim store: remote MCP + REST server (src/claimstore/mcp_http.py, the engine). Standard library only, no pip install.
# Build context is deploy/mcp/build/, staged by deploy/mcp/stage.sh (code + a snapshot of claim_store.sqlite).
FROM python:3.12-slim

WORKDIR /app
COPY claimstore/ /app/claimstore/
COPY config/ /app/config/
COPY claim_store.sqlite /app/claim_store.sqlite

ENV CLAIMSTORE_CONFIG=/app/config/project.toml \
    PDREP_STORE_DB=/app/claim_store.sqlite \
    PDREP_HTTP_HOST=0.0.0.0 \
    PORT=8080 \
    PYTHONUNBUFFERED=1

# uid 1000: Hugging Face Spaces run containers as this user
RUN useradd --uid 1000 --no-create-home app && chmod -R a-w /app
USER app
EXPOSE 8080
HEALTHCHECK --interval=30s --timeout=5s CMD python3 -c "import urllib.request,os;urllib.request.urlopen(f'http://127.0.0.1:{os.environ[\"PORT\"]}/health')"
CMD ["python3", "-I", "/app/claimstore/mcp_http.py"]

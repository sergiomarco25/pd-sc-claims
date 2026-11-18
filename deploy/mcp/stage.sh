#!/usr/bin/env bash
# Stage the container build context: the server code and a consistent snapshot of the claim store.
# Usage: deploy/mcp/stage.sh [path/to/claim_store.sqlite]
# Then:  docker build -t pd-claim-store deploy/mcp/build   (or podman / gcloud builds submit, see README.md)
set -euo pipefail
here="$(cd "$(dirname "$0")" && pwd)"
repo="$(cd "$here/../.." && pwd)"
db="${1:-${PDREP_STORE_DB:-$(python3 -I "$repo/src/claimstore/config.py" STORE_DB)}}"   # default from config/project.toml
out="$here/build"
rm -rf "$out"; mkdir -p "$out/claimstore" "$out/config"
cp "$here/Dockerfile" "$out/"
cp "$here/space_README.md" "$out/README.md"  # Hugging Face Space header (sdk: docker, app_port); harmless elsewhere
for f in __init__.py config.py kg.py mcp_server.py mcp_http.py; do cp "$repo/src/claimstore/$f" "$out/claimstore/"; done
# the project config: name, env prefix, server instructions and examples (paths in it are not used by the server)
cp "$repo/config/project.toml" "$repo/config/server_instructions.md" "$out/config/"
# sqlite backup API: a consistent copy even if the store is being rebuilt
/usr/bin/env python3 -I -c "import sqlite3,sys; s=sqlite3.connect(f'file:{sys.argv[1]}?mode=ro',uri=True); d=sqlite3.connect(sys.argv[2]); s.backup(d); d.close()" "$db" "$out/claim_store.sqlite"
echo "staged $out ($(du -sh "$out" | cut -f1)); store snapshot from $db"

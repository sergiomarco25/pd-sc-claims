"""Export the claim store as a knowledge graph (src/pdrep/kg.py).

Writes
  claims/store/kg_nodes.tsv, kg_edges.tsv           reviewable in git
  results/store/knowledge_graph.graphml             for Cytoscape / Gephi / networkx
  results/store/knowledge_graph.json                node-link JSON (networkx.node_link_graph compatible)

Standard library only: python3 -I scripts/57_knowledge_graph.py (any Python 3.9+). Run after 51b_claim_store.py.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep import kg  # noqa: E402


def main():
    g = kg.build()
    out = kg.DEFAULT_DB.parent
    g.write_tsv(ROOT / "claims" / "store" / "kg_nodes.tsv", ROOT / "claims" / "store" / "kg_edges.tsv")
    g.write_graphml(out / "knowledge_graph.graphml")
    (out / "knowledge_graph.json").write_text(json.dumps(g.to_json(), default=str))
    s = g.stats()
    print(f"{sum(s['nodes'].values())} nodes {s['nodes']}")
    print(f"{sum(s['edges'].values())} edges {s['edges']}")


if __name__ == "__main__":
    main()

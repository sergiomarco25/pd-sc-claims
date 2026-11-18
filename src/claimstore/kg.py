"""Knowledge graph view of the claim store (claims/STORE_DESIGN.md section 3, "graph export").

The store (results/store/claim_store.sqlite, built by scripts/51b_claim_store.py) stays the primary, relational record.
This module reads it and returns the same knowledge as typed nodes and edges, so that it can be drawn (atlas Graph
page), traversed (MCP server, scripts/mcp_server.py) or exported (scripts/57_knowledge_graph.py: TSV, GraphML, JSON).

Standard library only, so it runs with any Python 3.9+ (also outside the analysis environments).

Node ids are "<type>:<key>":
  claim:<claim_id>          published claim or finding of this study (attribute origin)
  paper:<paper_id>          source of claims; paper:this-study for our own analyses
  dataset:<dataset>         harmonised dataset (cohort)
  cell:<label>              cell type (CL id), cell:All = all cell types
  region:<label>            brain region (UBERON id)
  gene:<symbol>             gene named by a claim (signature members are not nodes; see signature_gene in the store)
  theme:<name>              biological theme
  statement:<CSxx>          consensus statement
  background:<BG-...>       reference knowledge from other evidence types (claims/background.yaml)
  question:<id>             question asked of the atlas (claims/questions.yaml)

Edges are directed (source, relation, target, attributes); traversal ignores direction.
"""
import heapq
import json
import math
import sqlite3
from collections import Counter, defaultdict
from pathlib import Path

from claimstore import config as cfg  # noqa: E402

DEFAULT_DB = cfg.STORE_DB

NODE_TYPES = ["claim", "paper", "dataset", "cell", "region", "gene", "theme", "statement", "background", "question"]
RELATIONS = {  # relation: (source type, target type, meaning)
    "asserts": ("paper", "claim", "the paper (or this study) makes the claim"),
    "about_cell_type": ("claim", "cell", "the claim is about this cell type"),
    "in_region": ("claim", "region", "the claim was observed in this region"),
    "names_gene": ("claim", "gene", "the claim names this gene"),
    "has_theme": ("claim", "theme", "the claim belongs to this biological theme"),
    "derived_from": ("claim", "dataset", "the claim was made from this dataset"),
    "tested_in": ("claim", "dataset", "the claim was re-tested in this dataset (attribute outcome)"),
    "supports_statement": ("claim", "statement", "the claim supports the consensus statement"),
    "opposes_statement": ("claim", "statement", "the claim opposes the consensus statement"),
    "agrees_with": ("claim", "claim", "two claims assert the same change (computed)"),
    "establishes": ("background", "statement", "reference knowledge establishes the statement"),
    "bears_on": ("background", "statement", "reference knowledge is relevant to the statement"),
    "supports": ("background", "claim", "reference knowledge supports the claim"),
    "consistent_with": ("background", "claim", "reference knowledge is consistent with the claim"),
    "contrasts_with": ("background", "claim", "reference knowledge contrasts with the claim"),
    "concerns_region": ("background", "region", "reference knowledge concerns this region"),
    "generated": ("paper", "dataset", "the paper generated the dataset"),
    "covers_region": ("dataset", "region", "the dataset sampled this region"),
    "cites": ("question", "claim", "the answer to the question cites the record"),
    "produced": ("question", "claim", "the question produced this finding"),
}


def _split(s):
    return [x.strip() for x in (s or "").split(";") if x and x.strip()]


def _json(s):
    """List stored as JSON text; falls back to ';'-separated text."""
    try:
        v = json.loads(s) if s else []
        return v if isinstance(v, list) else [v]
    except (TypeError, ValueError):
        return _split(s)


class Graph:
    def __init__(self):
        self.nodes = {}  # id -> dict(type, label, **attrs)
        self.edges = []  # (src, rel, dst, attrs)
        self.adj = defaultdict(list)  # id -> [(edge index)]
        self._seen = set()

    # ---------------------------------------------------------------- build
    def node(self, nid, type_, label, **attrs):
        if nid not in self.nodes:
            self.nodes[nid] = dict(type=type_, label=label, **{k: v for k, v in attrs.items() if v not in (None, "")})
        return nid

    def edge(self, src, rel, dst, **attrs):
        if src not in self.nodes or dst not in self.nodes or src == dst:
            return
        key = (src, rel, dst, attrs.get("role"))
        if key in self._seen:
            return
        self._seen.add(key)
        self.edges.append((src, rel, dst, {k: v for k, v in attrs.items() if v not in (None, "")}))
        i = len(self.edges) - 1
        self.adj[src].append(i)
        self.adj[dst].append(i)

    # ---------------------------------------------------------------- query
    def degree(self, nid):
        return len(self.adj.get(nid, ()))

    def resolve(self, name, types=None):
        """Node ids matching a free name: exact id, then exact key or label (any type), then substring of label."""
        if not name:
            return []
        name = str(name).strip()
        ok = (lambda nid: self.nodes[nid]["type"] in types) if types else (lambda nid: True)
        if name in self.nodes and ok(name):
            return [name]
        low = name.lower()
        exact = [n for n in self.nodes if ok(n) and (n.split(":", 1)[1].lower() == low or
                                                     str(self.nodes[n]["label"]).lower() == low)]
        if exact:
            return sorted(exact, key=lambda n: NODE_TYPES.index(self.nodes[n]["type"]))
        part = [n for n in self.nodes if ok(n) and (low in n.lower() or low in str(self.nodes[n]["label"]).lower())]
        return sorted(part, key=lambda n: (-self.degree(n), n))[:25]

    def incident(self, nid, relations=None):
        for i in self.adj.get(nid, ()):
            s, r, d, a = self.edges[i]
            if relations and r not in relations:
                continue
            yield (d if s == nid else s), r, ("out" if s == nid else "in"), a

    def neighbourhood(self, nid, hops=1, relations=None, types=None, limit=200):
        """BFS ball around nid: nodes (with hop distance) and the edges among them. Capped at `limit` nodes; nodes of
        higher degree are kept first within a hop when the cap is hit."""
        dist, frontier = {nid: 0}, [nid]
        for h in range(1, hops + 1):
            nxt = []
            for u in frontier:
                for v, r, _, _ in self.incident(u, relations):
                    if v not in dist and (not types or self.nodes[v]["type"] in types):
                        nxt.append(v)
            nxt = sorted(set(nxt), key=lambda v: -self.degree(v))
            room = limit - len(dist)
            for v in nxt[:max(room, 0)]:
                dist[v] = h
            frontier = nxt[:max(room, 0)]
            if room <= len(nxt):
                break
        keep = set(dist)
        edges = [e for e in self.edges if e[0] in keep and e[2] in keep and (not relations or e[1] in relations)]
        return dist, edges

    def _cheapest(self, a, b, max_hops, avoid_types, banned):
        """Dijkstra a -> b (undirected) where entering an intermediate node costs 1 + log10(degree), so specific links
        (a shared gene, a consensus statement) are preferred over hubs (all cell types, a large theme)."""
        best, heap = {a: 0.0}, [(0.0, 0, a, [])]
        while heap:
            cost, hops, u, path = heapq.heappop(heap)
            if u == b:
                return path
            if cost > best.get(u, float("inf")) or hops >= max_hops:
                continue
            if u != a and (self.nodes[u]["type"] in avoid_types or u in banned):
                continue
            for i in self.adj.get(u, ()):
                s, _, d, _ = self.edges[i]
                v = d if s == u else s
                c = cost + 1 + (0 if v == b else math.log10(max(self.degree(v), 1)))
                if c < best.get(v, float("inf")):
                    best[v] = c
                    heapq.heappush(heap, (c, hops + 1, v, path + [i]))
        return None

    def paths(self, a, b, max_hops=5, avoid_types=("theme",), k=3):
        """Up to k distinct cheapest paths a -> b. Each further path bans the intermediate nodes already used."""
        if a not in self.nodes or b not in self.nodes:
            return []
        out, banned = [], set()
        for _ in range(k):
            p = self._cheapest(a, b, max_hops, avoid_types, banned)
            if not p:
                break
            out.append([self.edges[i] for i in p])
            mids = {n for i in p for n in (self.edges[i][0], self.edges[i][2])} - {a, b}
            if not mids:
                break
            banned |= mids
        return out

    def stats(self):
        return dict(nodes=dict(Counter(n["type"] for n in self.nodes.values())),
                    edges=dict(Counter(e[1] for e in self.edges)))

    # ---------------------------------------------------------------- export
    def to_json(self):
        return dict(directed=True, multigraph=True,
                    nodes=[dict(id=k, **v) for k, v in self.nodes.items()],
                    links=[dict(source=s, relation=r, target=d, **a) for s, r, d, a in self.edges])

    def compact(self, similar=None, topics=None):
        """Small form for embedding in the atlas page: nodes [id, type, label, attrs], edges [i, j, rel, attrs] and,
        when given (similarity_links), similarity pairs [i, j, cosine] for the 3D 'by meaning' layout."""
        ids = list(self.nodes)
        ix = {n: i for i, n in enumerate(ids)}
        keep = ("verdict", "certainty", "origin", "paper", "n_donors", "short", "ontology_id", "n_claims", "region_group")
        nodes = [[n, self.nodes[n]["type"], str(self.nodes[n]["label"])[:220],
                  {k: self.nodes[n][k] for k in keep if k in self.nodes[n]}] for n in ids]
        edges = [[ix[s], ix[d], r, {k: a[k] for k in ("outcome", "role") if k in a}] for s, r, d, a in self.edges]
        sim = [[ix[a], ix[b], w] for a, b, w in (similar or []) if a in ix and b in ix]
        top = dict(labels=(topics or {}).get("labels", []),
                   members=[[ix[n], k] for n, k in ((topics or {}).get("members") or {}).items() if n in ix])
        return dict(nodes=nodes, edges=edges, similar=sim, topics=top,
                    relations={r: m for r, (_, _, m) in RELATIONS.items()})

    def write_tsv(self, nodes_path, edges_path):
        attr_keys = sorted({k for v in self.nodes.values() for k in v} - {"type", "label"})
        with open(nodes_path, "w") as f:
            f.write("\t".join(["id", "type", "label", "degree"] + attr_keys) + "\n")
            for n, v in sorted(self.nodes.items()):
                row = [n, v["type"], v["label"], str(self.degree(n))] + [str(v.get(k, "")) for k in attr_keys]
                f.write("\t".join(x.replace("\t", " ").replace("\n", " ") for x in row) + "\n")
        ek = sorted({k for *_, a in self.edges for k in a})
        with open(edges_path, "w") as f:
            f.write("\t".join(["source", "relation", "target"] + ek) + "\n")
            for s, r, d, a in sorted(self.edges, key=lambda e: (e[1], e[0], e[2])):
                f.write("\t".join([s, r, d] + [str(a.get(k, "")) for k in ek]) + "\n")

    def write_graphml(self, path):
        from xml.sax.saxutils import escape, quoteattr
        nk = sorted({k for v in self.nodes.values() for k in v})
        ek = ["relation"] + sorted({k for *_, a in self.edges for k in a})
        lines = ['<?xml version="1.0" encoding="UTF-8"?>',
                 '<graphml xmlns="http://graphml.graphdrawing.org/xmlns">']
        lines += [f'<key id="n_{k}" for="node" attr.name="{k}" attr.type="string"/>' for k in nk]
        lines += [f'<key id="e_{k}" for="edge" attr.name="{k}" attr.type="string"/>' for k in ek]
        lines.append('<graph id="pd_claims" edgedefault="directed">')
        for n, v in self.nodes.items():
            data = "".join(f'<data key="n_{k}">{escape(str(v[k]))}</data>' for k in nk if k in v)
            lines.append(f"<node id={quoteattr(n)}>{data}</node>")
        for s, r, d, a in self.edges:
            data = f'<data key="e_relation">{escape(r)}</data>' + "".join(
                f'<data key="e_{k}">{escape(str(a[k]))}</data>' for k in ek[1:] if k in a)
            lines.append(f"<edge source={quoteattr(s)} target={quoteattr(d)}>{data}</edge>")
        lines += ["</graph>", "</graphml>"]
        Path(path).write_text("\n".join(lines) + "\n")


EMBED_MODEL = "NeuML/pubmedbert-base-embeddings"  # as in scripts/56_claim_map.py


def similarity_links(db=None, k=6, min_sim=0.0, cache=None):
    """Nearest neighbours by meaning among claims, consensus statements and reference knowledge, from the sentence
    embeddings cached by scripts/56_claim_map.py (same texts, same keys). Returns [(node_a, node_b, cosine)], each pair
    once. Needs numpy (only here); records without a cached vector are skipped. Used for the 'by meaning' 3D layout.
    `cache`: the embedding file (default: claim_embeddings.npz next to the store; a release store reuses the full one,
    whose vectors are looked up by text)."""
    import hashlib
    import numpy as np
    cache = Path(cache) if cache else Path(db or DEFAULT_DB).parent / "claim_embeddings.npz"
    if not cache.exists():
        return []
    have = dict(np.load(cache, allow_pickle=True)["d"].item())
    con = connect(db)
    items = []
    for cid, st, ct, cs, rg in con.execute("SELECT claim_id, statement, cell_types, cell_state, regions FROM claim"):
        text = f"{st} Cell type: {(ct or '').replace(';', ', ')}." + (f" State: {cs}." if cs else "") + (
            f" Region: {rg.replace(';', ', ')}." if rg else "")
        items.append((f"claim:{cid}", text))
    items += [(f"statement:{i}", t) for i, t in con.execute("SELECT statement_id, statement FROM consensus_statement")]
    items += [(f"background:{i}", t) for i, t in con.execute("SELECT id, statement FROM background")]
    con.close()
    keyed = [(nid, have.get(hashlib.sha1(f"{EMBED_MODEL}|{t}".encode()).hexdigest())) for nid, t in items]
    keyed = [(nid, v) for nid, v in keyed if v is not None]
    if len(keyed) < 2:
        return []
    ids = [nid for nid, _ in keyed]
    E = np.vstack([v for _, v in keyed]).astype(np.float32)
    E /= np.linalg.norm(E, axis=1, keepdims=True)
    S = E @ E.T
    np.fill_diagonal(S, -1)
    out = {}
    for i in range(len(ids)):
        for j in np.argsort(-S[i])[:k]:
            if S[i, j] >= min_sim:
                a, b = sorted((i, int(j)))
                out[(a, b)] = float(S[i, j])
    return [(ids[a], ids[b], round(w, 3)) for (a, b), w in out.items()]


def topic_clusters(db=None, map_dir=None):
    """Topic clusters of claims (and statements / reference records) from the claim map (scripts/56_claim_map.py:
    claims/store/claim_map.tsv, claim_clusters.tsv), named by their dominant cell type and theme:
    "<cell type>: <theme>" when one cell type holds at least half the claims, else "<theme> (mixed cell types)";
    repeated names get the cluster's top TF-IDF terms. Returns dict(labels=[name per cluster], members={node_id: k})."""
    import csv
    d = Path(map_dir or Path(__file__).resolve().parents[2] / "claims" / "store")
    if not (d / "claim_map.tsv").exists():
        return dict(labels=[], members={})
    pref = {"claim": "claim", "statement": "statement", "reference": "background"}
    members = {}
    for r in csv.DictReader(open(d / "claim_map.tsv"), delimiter="\t"):
        if r.get("cluster") not in (None, "") and r["kind"] in pref:
            members[f"{pref[r['kind']]}:{r['id']}"] = int(float(r["cluster"]))
    terms = {int(r["cluster"]): r["label"] for r in csv.DictReader(open(d / "claim_clusters.tsv"), delimiter="\t")} \
        if (d / "claim_clusters.tsv").exists() else {}
    con = connect(db)
    cells, themes, size = defaultdict(Counter), defaultdict(Counter), Counter()
    for cid, cts, ths in con.execute("SELECT claim_id, cell_types, themes FROM claim"):
        k = members.get(f"claim:{cid}")
        if k is None:
            continue
        size[k] += 1
        cells[k].update(x for x in _split(cts) if x != "All")
        themes[k].update(_split(ths))
    con.close()
    n = max(members.values(), default=-1) + 1
    labels, seen = [], Counter()
    for k in range(n):
        th = themes[k].most_common(1)[0][0] if themes[k] else "unthemed"
        top = cells[k].most_common(1)
        name = f"{top[0][0]}: {th}" if top and top[0][1] >= size[k] / 2 else f"{th} (mixed cell types)"
        labels.append(name)
        seen[name] += 1
    for k, name in enumerate(list(labels)):  # repeated name: the largest cluster keeps it, others add distinctive terms
        if seen[name] < 2 or size[k] == max(size[j] for j, x in enumerate(labels) if x == name) or not terms.get(k):
            continue
        cand = sorted([t for t in terms[k].split(" · ") if t.lower() not in name.lower()], key=lambda t: -len(t.split()))
        picked = []
        for t in cand:
            if not any(t in q for q in picked) and len(picked) < 2:
                picked.append(t)
        labels[k] = f"{name} · {', '.join(picked)}" if picked else name
    return dict(labels=labels, members=members)


def connect(db=None):
    return sqlite3.connect(f"file:{db or DEFAULT_DB}?mode=ro", uri=True)


def build(db=None):
    """Read the claim store and return its knowledge graph."""
    con = connect(db)
    con.row_factory = sqlite3.Row
    q = lambda sql: [dict(r) for r in con.execute(sql)]  # noqa: E731
    g = Graph()

    vocab_r = {r["label"]: r for r in q("SELECT * FROM vocab_regions")}
    vocab_c = {r["label"]: r for r in q("SELECT * FROM vocab_cell_types")}
    src = {r["paper_id"]: r for r in q("SELECT * FROM source")}

    def paper(pid):
        if pid == "this-study":
            return g.node("paper:this-study", "paper", "This study (our analyses)", kind="our_analysis")
        s = src.get(pid, {})
        lab = f"{s.get('first_author') or pid} {s.get('year') or ''}".strip()
        return g.node(f"paper:{pid}", "paper", lab, title=s.get("title"), doi=s.get("doi"),
                      journal=s.get("journal_or_preprint"), kind=s.get("kind") or "paper")

    def region_label(text):
        """Free region text (e.g. 'vagal nucleus (DMNX)', 'motor cortex') -> vocabulary label, or the text itself."""
        t = text.strip()
        for lab, v in vocab_r.items():
            names = {lab.lower(), str(v.get("ontology_label") or "").lower(), f"{lab} cortex".lower()}
            if t.lower() in names or t.endswith(f"({lab})"):
                return lab
        return t

    def region(lab):
        v = vocab_r.get(lab, {})
        return g.node(f"region:{lab}", "region", lab, ontology_id=v.get("ontology_id"),
                      region_group=v.get("region_group"))

    def cell(lab):
        v = vocab_c.get(lab, {})
        return g.node(f"cell:{lab}", "cell", "All cell types" if lab == "All" else lab,
                      ontology_id=v.get("ontology_id"))

    # datasets
    for d in q("SELECT * FROM dataset"):
        g.node(f"dataset:{d['dataset']}", "dataset", d["cohort"] if d["cohort"] != d["dataset"] else d["dataset"],
               modality=d["modality"], n_donors=d["n_donors"], regions=d["regions"], accession=d["accession"],
               donors_by_disease=d["donors_by_disease"], ancestry=d["ancestry_summary"])
        g.edge(paper(d["paper"]), "generated", f"dataset:{d['dataset']}")
        for r in _split(d["regions"]):
            g.edge(f"dataset:{d['dataset']}", "covers_region", region(r))

    # claims
    for c in q("SELECT * FROM claim"):
        cid = f"claim:{c['claim_id']}"
        g.node(cid, "claim", c["statement"], paper=c["paper"], origin=c["origin"], verdict=c["verdict"] or "Untested",
               certainty=c["certainty"], rating=c["rating"], claim_type=c["claim_type"], direction=c["direction"],
               predicate=c["predicate"], subject=c["subject"], contrast=c["contrast"], centrality=c["centrality"],
               cell_state=c["cell_state"], tier=c["tier"])
        g.edge(paper(c["paper"]), "asserts", cid)
        for x in _split(c["cell_types"]):
            g.edge(cid, "about_cell_type", cell(x))
        for x in _split(c["regions"]):
            g.edge(cid, "in_region", region(x))
        for x in _split(c["genes"]):
            g.edge(cid, "names_gene", g.node(f"gene:{x}", "gene", x))
        for x in _split(c["themes"]):
            g.edge(cid, "has_theme", g.node(f"theme:{x}", "theme", x))
        for x in _split(c["source_dataset"]):
            g.edge(cid, "derived_from", f"dataset:{x}")

    for e in q("SELECT claim_id, dataset, role, outcome, effect, z, p, region_match FROM claim_dataset "
               "WHERE outcome IS NOT NULL AND outcome != 'not_testable'"):
        g.edge(f"claim:{e['claim_id']}", "tested_in", f"dataset:{e['dataset']}", role=e["role"], outcome=e["outcome"],
               effect=e["effect"], z=e["z"], p=e["p"], region_match=e["region_match"])

    for s in q("SELECT * FROM consensus_statement"):
        g.node(f"statement:{s['statement_id']}", "statement", s["statement"], short=s["short"], verdict=s["verdict"])
    for r in q("SELECT * FROM claim_relation"):
        b = f"statement:{r['claim_b']}" if r["relation"].endswith("_statement") else f"claim:{r['claim_b']}"
        g.edge(f"claim:{r['claim_a']}", r["relation"], b, basis=r["basis"], note=r["note"])

    for b in q("SELECT * FROM background"):
        bid = g.node(f"background:{b['id']}", "background", b["statement"], evidence_type=b["evidence_type"],
                     measurement_note=b["measurement_note"])
        est = set(_json(b["establishes"]))
        for s in _json(b["statements"]):
            g.edge(bid, "establishes" if s in est else "bears_on", f"statement:{s}")
        for rel, col in (("supports", "supports_claims"), ("consistent_with", "consistent_with"),
                         ("contrasts_with", "contrasts_with")):
            for x in _json(b[col]):
                g.edge(bid, rel, f"claim:{x}")
        for r in _json(b["regions"]):
            g.edge(bid, "concerns_region", region(region_label(r)))

    for qu in q("SELECT * FROM question"):
        qid = g.node(f"question:{qu['id']}", "question", qu["text"], status=qu["status"], answer=qu["answer"],
                     asked_by=qu["asked_by"], date=qu["date"])
        for rec in _json(qu["records"]):
            kind, _, key = str(rec).partition(":")
            tgt = {"claim": f"claim:{key}", "stmt": f"statement:{key}"}.get(kind)
            if tgt:
                g.edge(qid, "cites", tgt)
        prop = qu["proposed"]
        try:
            fid = (json.loads(prop) if prop else {}).get("finding", {}).get("id")
        except (TypeError, ValueError, AttributeError):
            fid = None
        if fid:
            g.edge(qid, "produced", f"claim:{fid}")

    for n, v in g.nodes.items():  # claim counts on hub nodes, for sizing
        if v["type"] != "claim":
            v["n_claims"] = sum(1 for u, *_ in g.incident(n) if g.nodes[u]["type"] == "claim")
    info = q("SELECT * FROM build_info")
    g.info = info[0] if info else {}
    con.close()
    return g

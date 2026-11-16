"""Claim map: every claim, consensus statement and reference-knowledge record placed in one 2D space by meaning.

Text per item = statement + cell types + cell state + regions; embedded with a biomedical sentence model
(NeuML/pubmedbert-base-embeddings, cached in results/store/claim_embeddings.npz by text), laid out with UMAP
(cosine) and clustered with Leiden on a k-nearest-neighbour graph of the claims. Clusters are named by their most
distinctive words (class-based TF-IDF). Statements and reference records are placed in the same space but do not
define clusters.

Outputs: claims/store/claim_map.tsv (id, kind, x, y, cluster), claims/store/claim_clusters.tsv (cluster, label, n,
x, y, top themes). The atlas map page reads both (51c).

Run (needs sentence-transformers, installed under scratch/pylib):
  PYTHONPATH=scratch/pylib HF_HOME=scratch/hf python scripts/56_claim_map.py [--db STORE --out DIR]
--db maps another store (e.g. a public release) and --out writes its map elsewhere (the public atlas uses
results/public/claim_map); embeddings always come from and go to the cache of the configured store.
"""
import hashlib
import re
import sqlite3
from pathlib import Path

import igraph as ig
import leidenalg
import numpy as np
import pandas as pd
import umap
import yaml
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.neighbors import NearestNeighbors

ROOT = Path(__file__).resolve().parents[1]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

STORE = cfg.RESULTS / "store"
MODEL = "NeuML/pubmedbert-base-embeddings"
STOP = set("pd parkinson parkinson's disease patients donors control controls compared higher lower increased decreased "
           "reduced expression cells cell type types brain human pd's vs versus across within show shows".split())


def items(db=None):
    con = sqlite3.connect(db or STORE / "claim_store.sqlite")
    c = pd.read_sql("SELECT claim_id id, statement, cell_types, cell_state, regions, themes FROM claim", con)
    s = pd.read_sql("SELECT statement_id id, statement FROM consensus_statement", con)
    con.close()
    c["kind"] = "claim"
    c["text"] = c.apply(lambda r: f"{r.statement} Cell type: {r.cell_types.replace(';', ', ')}."
                        + (f" State: {r.cell_state}." if r.cell_state else "")
                        + (f" Region: {r.regions.replace(';', ', ')}." if r.regions else ""), axis=1)
    s["kind"], s["text"], s["themes"] = "statement", s.statement, ""
    bg = yaml.safe_load((ROOT / "claims" / "background.yaml").read_text()) or []
    b = pd.DataFrame([{"id": x["id"], "statement": x["statement"], "text": x["statement"], "kind": "reference", "themes": ""}
                      for x in bg])
    return pd.concat([c, s, b], ignore_index=True).fillna("")


def embed(texts):
    cache = STORE / "claim_embeddings.npz"
    keys = [hashlib.sha1(f"{MODEL}|{t}".encode()).hexdigest() for t in texts]
    have = dict(np.load(cache, allow_pickle=True)["d"].item()) if cache.exists() else {}
    todo = [i for i, k in enumerate(keys) if k not in have]
    if todo:
        from sentence_transformers import SentenceTransformer
        m = SentenceTransformer(MODEL)
        vecs = m.encode([texts[i] for i in todo], normalize_embeddings=True, batch_size=32, show_progress_bar=False)
        for i, v in zip(todo, vecs):
            have[keys[i]] = v.astype(np.float32)
        np.savez_compressed(cache, d=np.array(have, dtype=object))
    print(f"embedded {len(todo)} new, {len(texts) - len(todo)} from cache")
    return np.vstack([have[k] for k in keys])


def label_clusters(df):
    docs = df.groupby("cluster").text.apply(" ".join)
    cv = CountVectorizer(token_pattern=r"[A-Za-z][A-Za-z0-9\-]{2,}", lowercase=True, ngram_range=(1, 2),
                         stop_words=list(STOP | set(CountVectorizer(stop_words="english").get_stop_words())))
    X = cv.fit_transform(docs).toarray().astype(float)
    tf = X / X.sum(1, keepdims=True)
    idf = np.log(1 + X.sum(1).mean() / (X.sum(0) + 1))
    score = tf * idf
    vocab = np.array(cv.get_feature_names_out())
    out = {}
    for i, cl in enumerate(docs.index):
        top = [w for w in vocab[np.argsort(-score[i])] if not any(w in o or o in w for o in [])][:12]
        picked = []
        for w in top:  # drop words already covered by a picked bigram, keep three
            if any(w in p for p in picked):
                continue
            picked.append(w)
            if len(picked) == 3:
                break
        out[cl] = " · ".join(picked)
    return out


def main():
    import argparse
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--db", help="store to map (default: the configured store)")
    ap.add_argument("--out", help="folder for claim_map.tsv and claim_clusters.tsv (default: claims/store)")
    a = ap.parse_args()
    df = items(a.db)
    E = embed(df.text.tolist())
    xy = umap.UMAP(n_neighbors=15, min_dist=0.35, metric="cosine", random_state=0).fit_transform(E)
    df["x"], df["y"] = xy[:, 0].round(4), xy[:, 1].round(4)
    cl_mask = df.kind.eq("claim").values
    Ec = E[cl_mask]
    k = 12
    nn = NearestNeighbors(n_neighbors=k + 1, metric="cosine").fit(Ec)
    dist, idx = nn.kneighbors(Ec)
    edges, w = [], []
    for i in range(len(Ec)):
        for j, d in zip(idx[i, 1:], dist[i, 1:]):
            edges.append((i, int(j))); w.append(1 - d)
    g = ig.Graph(n=len(Ec), edges=edges, directed=False)
    g.es["weight"] = w
    g.simplify(combine_edges="max")
    part = leidenalg.find_partition(g, leidenalg.RBConfigurationVertexPartition, weights="weight",
                                    resolution_parameter=1.0, seed=0)
    df["cluster"] = -1
    df.loc[cl_mask, "cluster"] = np.array(part.membership)
    # statements and reference records take the cluster of their nearest claims
    others = ~cl_mask
    if others.any():
        nn2 = NearestNeighbors(n_neighbors=5, metric="cosine").fit(Ec)
        _, idx2 = nn2.kneighbors(E[others])
        memb = np.array(part.membership)
        df.loc[others, "cluster"] = [int(pd.Series(memb[r]).mode()[0]) for r in idx2]
    labels = label_clusters(df[df.kind.eq("claim")])
    cl = (df[df.kind.eq("claim")].groupby("cluster")
          .agg(n=("id", "size"), x=("x", "mean"), y=("y", "median"),
               themes=("themes", lambda s: "; ".join(pd.Series([t for v in s for t in str(v).split(";") if t])
                                                   .value_counts().head(2).index)))
          .reset_index())
    cl["label"] = cl.cluster.map(labels)
    out = Path(a.out) if a.out else ROOT / "claims" / "store"
    out.mkdir(parents=True, exist_ok=True)
    df[["id", "kind", "x", "y", "cluster"]].to_csv(out / "claim_map.tsv", sep="\t", index=False)
    cl[["cluster", "label", "n", "x", "y", "themes"]].round(3).to_csv(out / "claim_clusters.tsv", sep="\t", index=False)
    print(f"{cl_mask.sum()} claims, {others.sum()} statements/references, {len(cl)} clusters")
    print(cl.sort_values("n", ascending=False)[["cluster", "n", "label", "themes"]].to_string(index=False))


if __name__ == "__main__":
    main()

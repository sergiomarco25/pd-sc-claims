"""The atlas: one static HTML page (works offline) for browsing a claim store: claims with their ratings and re-tests,
consensus statements, reference knowledge, datasets, genes, a knowledge graph and an "Ask" box.

Usage: python3 -I src/claimstore/cli.py atlas [--public] [--db STORE] [--out FILE]

Reads the store (results/store/claim_store.sqlite, or --db: the claim_store.sqlite of a public release) and writes
reports/atlas.html (--public: reports/atlas_public.html, with questions pointed at the public repository's issue form)
plus <out>.artifact.html, the same page without the html/head skeleton, for publishing as a claude.ai Artifact (with
the `sample` capability, "Ask the atlas" can call Claude).

config/project.toml [atlas], all optional:
  title, headline, lead, crumb   page title and overview wording (default: [project] title and description)
  disease, case_label            names for the case group in captions (its short name; "<name> vs control")
  scope                          what the records are about, for the "Ask" prompt
  examples, ask_placeholder      example questions; stop_words for search
  region_groups, region_group_of columns of the data landscape (default: one per region)
  key_findings                   [{label, title, text, href}] cards on the overview
  graph_start                    knowledge-graph node ids to start from
  template                       a project's own page (path from the repository); default the engine's
  extras                         "module:function(data, public)" returning tables a project adds (composition,
                                 programmes, regions, specificity, genetics, genes, pathways; see the template)
  out                            output path (default reports/atlas.html)
The engine's page shows the biology sections only when a project's extras provide their tables.
"""
import argparse
import html as _html
import importlib
import json
import re
import sqlite3
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from claimstore import config as cfg
from claimstore import kg

ROOT = cfg.REPO
ATLAS = cfg.CONFIG.get("atlas", {})
TEMPLATE = Path(__file__).with_name("atlas_template.html")
FIELD_KEYS = ("title", "headline", "lead", "crumb", "disease", "case_label", "scope", "examples", "ask_placeholder",
              "stop_words", "region_groups", "region_group_of", "key_findings", "graph_start", "regions_label")


def recs(df, digits=4):
    df = df.copy()
    for c in df.select_dtypes("number"):
        df[c] = df[c].round(digits)
    return df.replace({np.nan: None}).to_dict("records")


def read_store(db):
    """Claims, datasets, tests, relations, consensus statements, papers, build info, gene -> claims and reference
    knowledge from a claim store."""
    con = sqlite3.connect(db)
    q = lambda s: pd.read_sql(s, con)  # noqa: E731
    claims = q("""SELECT claim_id id, paper, source_dataset ds, statement, quote, location, organisms, organism_ids,
                  organism_note, sample_type, region_asserted, regions, region_ids, region_groups, cell_types,
                  cell_type_ids, cell_state, assay, assay_id, subject_kind, subject, predicate, direction, contrast,
                  contrast_kind, case_group, reference_group, comparison_raw, claim_type, centrality,
                  unit_of_analysis, reported_stat, method, genes, needs_review, context_text, evidence_text, tier,
                  rank, themes, verdict, certainty, rating, family, start, k_cohorts, est, lo, hi, i2, concordance,
                  z_comb, detect, donors, rob, inconsistency, indirectness, imprecision, dependence, upgrade, reasons,
                  not_assessed, scope_tested, scope_held, caveats, origin, analysis_script, analysis_results,
                  analysis_date, notes FROM claim""")
    claims["rank"] = pd.to_numeric(claims["rank"], errors="coerce")
    ds = q("SELECT * FROM dataset")
    ev = q("""SELECT claim_id c, dataset d, role, variant v, region_match rm, outcome o, effect, z, p,
              n_case nc, n_ref nr, regions_tested reg, classes_tested cls, note FROM claim_dataset""")
    for k in ("effect", "z", "p"):
        ev[k] = pd.to_numeric(ev[k], errors="coerce")
    rel = q("SELECT claim_a a, relation r, claim_b b, basis, note FROM claim_relation")
    cs = q("SELECT * FROM consensus_statement")
    bgt = q("SELECT * FROM background")
    for col in ("references", "statements", "establishes", "supports_claims", "consistent_with", "contrasts_with", "regions"):
        bgt[col] = bgt[col].map(lambda v: json.loads(v) if v else [])
    background = bgt.replace({np.nan: None}).to_dict("records")
    sig = q("""SELECT cs.claim_id c, sg.gene g FROM claim_signature cs
               JOIN signature_gene sg ON cs.signature_id = sg.signature_id""").dropna().drop_duplicates()
    papers = q("SELECT paper_id id, first_author, year, title, doi, journal_or_preprint journal FROM source "
               "WHERE paper_id IN (SELECT DISTINCT paper FROM claim)")
    info = q("SELECT * FROM build_info").iloc[0].to_dict()
    src = q("SELECT paper_id, first_author, year FROM source").set_index("paper_id")
    con.close()
    named = {}
    for _, r in claims.iterrows():
        for g in filter(None, (r.genes or "").split(";")):
            named.setdefault(g, set()).add(r.id)
    in_sig = sig.groupby("g").c.apply(set).to_dict()
    gene_claims = {g: [sorted(named.get(g, set())), sorted(in_sig.get(g, set()) - named.get(g, set()))]
                   for g in set(named) | set(in_sig)}
    ds["label"] = ds.apply(lambda d: dataset_label(d, src), axis=1) if len(ds) else []
    return dict(claims=claims, datasets=ds, evidence=ev, relations=rel, statements=cs, papers=papers, info=info,
                gene_claims=gene_claims, background=background)


def dataset_label(d, src):
    """'First author year', with the cohort when it differs from the dataset."""
    a = (f"{src.first_author.get(d.dataset, '')} {int(float(src.year.get(d.dataset)))}"
         if d.dataset in src.index and src.year.get(d.dataset) else d.dataset)
    return a if d.cohort == d.dataset else f"{a} · {d.cohort}"


def stats(claims, ds, nuclei=None, nuclei_source=None):
    tested = claims.verdict.fillna("Untested").ne("Untested")
    held = claims.verdict.fillna("").str.startswith("Supported")
    strong = held & claims.certainty.isin(["High", "Moderate"])
    regions = sorted({r for x in ds.regions for r in str(x).split(";")} - {"other"})
    findings = claims.origin.eq("finding")
    return dict(findings=int(findings.sum()), papers=int(claims.paper[~findings].nunique()), datasets=len(ds),
                donors=int(pd.to_numeric(ds.n_donors, errors="coerce").fillna(0).sum()), regions=len(regions),
                claims=int((~findings).sum()), tested=int((tested & ~findings).sum()),
                held=int((held & ~findings).sum()), strong=int((strong & ~findings).sum()),
                nuclei=nuclei, nuclei_source=nuclei_source)


def question_form():
    """GitHub issue form for questions: of the public repository (config [release] public_repository), else of the
    repository in CITATION.cff."""
    if cfg.PUBLIC_REPOSITORY:
        return cfg.PUBLIC_REPOSITORY.rstrip("/") + "/issues/new?template=question.yml"
    f = ROOT / "CITATION.cff"
    m = re.search(r'^repository-code:\s*"?([^"\n]+)"?', f.read_text(), re.M) if f.exists() else None
    return m.group(1).rstrip("/") + "/issues/new?template=question.yml" if m else None


def knowledge_graph(db, public=False):
    """The knowledge graph in the page's compact form, with the 'by meaning' links and topic clusters when their
    inputs exist (claim embeddings next to the store; the claim map)."""
    return kg.build(db).compact(kg.similarity_links(db, cache=cfg.STORE_DB.parent / "claim_embeddings.npz"),
                                kg.topic_clusters(db, map_dir=cfg.RESULTS / "public" / "claim_map" if public else None))


def contrast_labels(db):
    """Contrast code -> label, from the store's contrast vocabulary."""
    try:
        con = sqlite3.connect(db)
        v = pd.read_sql("SELECT * FROM vocab_contrasts", con)
        con.close()
    except Exception:
        return {}
    code = next((c for c in ("code", "contrast", "id") if c in v), None)
    lab = next((c for c in ("label", "name") if c in v), None)
    return dict(zip(v[code], v[lab])) if code and lab else {}


def harmonise_classes():
    """The cell classes of config/harmonise.toml, in their order (none if the project has no such file)."""
    import tomllib
    f = ROOT / "config" / "harmonise.toml"
    return list(tomllib.loads(f.read_text()).get("cell_classes", [])) if f.exists() else []


def field():
    p = cfg.PROJECT
    f = {k: ATLAS[k] for k in FIELD_KEYS if k in ATLAS}
    f.setdefault("title", p.get("title", "Claim atlas"))
    f.setdefault("lead", p.get("description", ""))
    return f


def data(db, public=False):
    """Everything the engine's page shows, before a project's extras."""
    s = read_store(db)
    claims = s["claims"]
    no_theme = claims.themes.fillna("").eq("")  # projects without a theme table: the claim type groups claims
    claims.loc[no_theme, "themes"] = claims.claim_type[no_theme].fillna("other").str.replace("_", " ")
    cells = [c for c in harmonise_classes() if c != "Unknown"] or sorted(
        {c for x in claims.cell_types.dropna() for c in x.split(";") if c and c != "All"})
    return dict(
        claims=recs(claims), datasets=recs(s["datasets"]), evidence=recs(s["evidence"], 3),
        relations=recs(s["relations"]), statements=recs(s["statements"]), background=s["background"],
        papers=recs(s["papers"]), gene_claims=s["gene_claims"], info={k: str(v) for k, v in s["info"].items()},
        stats=stats(claims, s["datasets"]), cell_order=cells + ["ALL"], field=field(), contrast_labels=contrast_labels(db),
        public=bool(public), links={"question_form": question_form()} if public else {},
        kg=knowledge_graph(db, public))


def render(d, template=TEMPLATE):
    f = d.get("field", {})
    page = Path(template).read_text()
    page = page.replace("__TITLE__", _html.escape(f.get("title", "Claim atlas")))
    page = page.replace("__BRAND__", _html.escape(f.get("brand", f.get("title", "Claim atlas"))))
    return page.replace("/*__DATA__*/null", json.dumps(d, separators=(",", ":"), default=str))


def write(page, out, artifact=True):
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(page)
    if artifact:  # claude.ai Artifact version: the publisher adds the doctype/head/body skeleton
        head, body = page.split("<body>", 1)
        art = head[head.index("<title>"):head.index("</head>")] + body.rsplit("</body>", 1)[0]
        out.with_name(out.stem + ".artifact.html").write_text(art)
    return out


def _call(spec, *args):
    mod, fn = spec.split(":")
    return getattr(importlib.import_module(mod), fn)(*args)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--public", action="store_true", help="public atlas (questions go to the public issue form)")
    ap.add_argument("--db", help="claim store (default: the project's; for --public, a release's store)")
    ap.add_argument("--out", help="output HTML")
    a = ap.parse_args(sys.argv[1:] if argv is None else argv)
    db = Path(a.db) if a.db else cfg.STORE_DB
    if not db.exists():
        sys.exit(f"no claim store at {db}: run `cli.py store` first")
    d = data(db, a.public)
    if ATLAS.get("extras"):
        d.update(_call(ATLAS["extras"], d, a.public) or {})
    tpl = ROOT / ATLAS["template"] if ATLAS.get("template") else TEMPLATE
    default = ATLAS.get("out", "reports/atlas.html")
    out = Path(a.out) if a.out else ROOT / (default.replace(".html", "_public.html") if a.public else default)
    write(render(d, tpl), out, artifact=not a.public)
    print(f"{out} ({out.stat().st_size / 1e6:.1f} MB): {len(d['claims'])} claims, {len(d['datasets'])} datasets, "
          f"{len(d['evidence'])} tests, {len(d['kg']['nodes'])} graph nodes")
    return 0


if __name__ == "__main__":
    sys.exit(main())

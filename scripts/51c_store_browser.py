"""PD single-cell atlas: a visual, browsable summary of what changes in the PD brain, linked to claims and data.

Reads the claim store (results/store/claim_store.sqlite, 51b_claim_store.py) and the cross-dataset biology results
(results/biology/*: composition, programme, region, specificity, genetics, gene meta-analyses), embeds them as one
compact JSON in scripts/atlas_template.html and writes reports/pd_atlas.html (static, works offline) and
reports/pd_atlas.artifact.html (same page without the html/head skeleton, for publishing as a claude.ai Artifact with
capabilities {sample: {}} so "Ask the atlas" can call Claude).

Run: ~/miniconda3/envs/rapids_singlecell/bin/python scripts/51c_store_browser.py
Public atlas (claims/COMMUNITY_PLAN.md step 5.1): ... 51c_store_browser.py --public --db results/releases/v<V>/claim_store.sqlite
  reads the public release store and the biology tables recomputed without private datasets (results/public/biology,
  results/public/psychad), leaves out the private papers' gene modules ([release] private_modules), the manuscript's
  model figure and the nuclei count, points questions to the GitHub issue form, and writes reports/pd_atlas_public.html.
"""
import argparse
import base64
import json
import sqlite3
from io import BytesIO
from pathlib import Path

import numpy as np
import pandas as pd

import sys
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from pdrep import kg  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]

from pdrep import config as cfg  # noqa: E402

RES = cfg.RESULTS
DB, BIO, PSY = RES / "store" / "claim_store.sqlite", RES / "biology", RES / "psychad"

BG_RECORDS = []  # reference knowledge (claims/background.yaml via the store)
# programme -> (group, short display name)
PROGRAMMES = {
    "heat response": ("Core stress response", "Heat-shock response"),
    "chaperone folding": ("Core stress response", "Chaperone folding"),
    "unfolded protein response": ("Core stress response", "Unfolded protein response"),
    "hypoxia": ("Core stress response", "Hypoxia"),
    "cold-inducible (CIRBP, RBM3)": ("Core stress response", "Cold-inducible RNA-binding (CIRBP, RBM3)"),
    "oxidative phosphorylation": ("Metabolism", "Oxidative phosphorylation"),
    "UPR astrocytes (Barba-Reyes)": ("Published stress states", "UPR astrocytes (Barba-Reyes)"),
    "UPR microglia (Barba-Reyes)": ("Published stress states", "UPR microglia (Barba-Reyes)"),
    "UPR astrocytes, Astrocytes3 (Martirosyan)": ("Published stress states", "UPR astrocytes, Astrocytes3 (Martirosyan)"),
    "disease-associated oligodendrocytes (Adams)": ("Published stress states", "Disease-associated oligodendrocytes (Adams)"),
    "stress oligodendrocytes, Oligos5 (Martirosyan)": ("Published stress states", "Stress oligodendrocytes, Oligos5 (Martirosyan)"),
    "HSP90AA1-high microglia (Smajić)": ("Published stress states", "HSP90AA1-high microglia (Smajić)"),
    "GFAP/MT astrocytes, Astrocytes4 (Martirosyan)": ("Published stress states", "GFAP/MT astrocytes, Astrocytes4 (Martirosyan)"),
    "reactive astrocytes (Barba-Reyes)": ("Published inflammatory / activation states", "Reactive astrocytes (Barba-Reyes)"),
    "AR2 astrocytes (Barba-Reyes)": ("Published inflammatory / activation states", "AR2 astrocytes (Barba-Reyes)"),
    "reactive microglia (Barba-Reyes)": ("Published inflammatory / activation states", "Reactive microglia (Barba-Reyes)"),
    "DAM microglia (Barba-Reyes)": ("Published inflammatory / activation states", "DAM microglia (Barba-Reyes)"),
    "GPNMB-high microglia (Smajić)": ("Published inflammatory / activation states", "GPNMB-high microglia (Smajić)"),
    "activated microglia, Micr-1 (Xu)": ("Published inflammatory / activation states", "Activated microglia, Micr-1 (Xu)"),
    "activated astrocytes, Ast-1 (Xu)": ("Published inflammatory / activation states", "Activated astrocytes, Ast-1 (Xu)"),
    "C3/CD44 astrocytes, Astrocytes1 (Martirosyan)": ("Published inflammatory / activation states", "C3/CD44 astrocytes, Astrocytes1 (Martirosyan)"),
    "activated myeloid cells (Ma)": ("Published inflammatory / activation states", "Activated myeloid cells (Ma)"),
    "inflammatory response": ("Inflammation pathways", "Inflammatory response"),
    "TNFa / NF-kB": ("Inflammation pathways", "TNFα / NF-κB"),
    "interferon gamma": ("Inflammation pathways", "Interferon γ"),
    "SOX6/AGTR1 DA subtype (Kamath)": ("Neuronal subtype", "SOX6/AGTR1 DA subtype (Kamath)"),
}
CELL_ORDER = ["DA neurons", "Excitatory neurons", "Inhibitory neurons", "MSNs", "Neurons (other)", "Astrocytes",
              "Microglia", "Oligodendrocytes", "OPCs", "Endothelial", "Pericytes", "Fibroblasts", "T cells",
              "Immune (other)", "Ependymal", "ALL"]
# brain schematic (midsagittal, front left) from scripts/39_region_map.py; Braak Lewy stage at which a site is reached
SITES = {"prefrontal cortex": (1.6, 6.1, 5), "motor cortex": (4.8, 7.55, 6), "cingulate cortex": (4.4, 5.8, 4),
         "visual cortex": (8.6, 5.0, 6), "temporal cortex": (3.9, 2.75, 4), "amygdala": (2.8, 3.45, 4),
         "caudate": (3.3, 5.0, None), "putamen": (3.85, 4.3, None), "globus pallidus": (4.45, 3.95, None),
         "substantia nigra": (5.35, 3.25, 3), "vagal nucleus (DMNX)": (6.0, 1.15, 1), "cerebellum": (7.7, 2.0, None)}
# label offset (dx, dy in schematic units) and text anchor per site, as in 39_region_map.py
SITE_LAB = {"prefrontal cortex": (-0.1, 0.45, "middle"), "motor cortex": (0, 0.45, "middle"), "cingulate cortex": (0.4, 0.15, "start"),
            "visual cortex": (0, 0.4, "middle"), "temporal cortex": (0, -0.45, "middle"), "amygdala": (-0.35, 0, "end"),
            "caudate": (-0.35, 0.1, "end"), "putamen": (-0.35, 0.05, "end"), "globus pallidus": (0.35, 0.3, "start"),
            "substantia nigra": (0.4, 0, "start"), "vagal nucleus (DMNX)": (0.4, 0, "start"), "cerebellum": (0.35, 0.3, "start")}
SITE_REGIONS = {"prefrontal cortex": ["PFC", "DLPFC"], "motor cortex": ["motor"], "cingulate cortex": ["cingulate"],
                "visual cortex": ["visual"], "temporal cortex": ["temporal", "MTG"], "amygdala": ["amygdala"],
                "caudate": ["caudate"], "putamen": ["putamen"], "globus pallidus": ["GPi"],
                "substantia nigra": ["SN", "midbrain"], "vagal nucleus (DMNX)": ["DMNX"], "cerebellum": ["cerebellum"]}


ATLAS = ROOT / "scripts" / "atlas_assets" / "ggseg"
# site -> atlas parcels (ggseg Desikan-Killiany cortex, FreeSurfer aseg); "clip" bands split the brainstem by level
ANAT_MAP = {
    "prefrontal cortex": [("dk", "superiorfrontal_left"), ("dk", "rostralmiddlefrontal_left"), ("dk", "caudalmiddlefrontal_left")],
    "motor cortex": [("dk", "precentral_left"), ("dk", "paracentral_left")],
    "temporal cortex": [("dk", "superiortemporal_left"), ("dk", "middletemporal_left"), ("dk", "inferiortemporal_left")],
    "cingulate cortex": [("dk", "rostralanteriorcingulate_left"), ("dk", "caudalanteriorcingulate_left"), ("dk", "posteriorcingulate_left"),
                         ("dk", "isthmuscingulate_left")],
    "visual cortex": [("dk", "pericalcarine_left")],
    "caudate": [("aseg", "Left-Caudate"), ("aseg", "Right-Caudate")],
    "putamen": [("aseg", "Left-Putamen"), ("aseg", "Right-Putamen")],
    "globus pallidus": [("aseg", "Left-Pallidum"), ("aseg", "Right-Pallidum")],
    "amygdala": [("aseg", "Left-Amygdala"), ("aseg", "Right-Amygdala")],
    "cerebellum": [("aseg", "Cerebellum-Cortex"), ("aseg", "Cerebellum-White-Matter")],
    "substantia nigra": [("aseg", "Brain-Stem", (None, 790))],      # midbrain level of the brainstem (top of the outline)
    "vagal nucleus (DMNX)": [("aseg", "Brain-Stem", (1000, None))],  # medulla level (lowest part)
}
ANAT_LABEL = {"substantia nigra": "midbrain · substantia nigra", "vagal nucleus (DMNX)": "medulla · vagal nucleus (DMNX)", "visual cortex": "primary visual cortex"}
VIEWS = [("lateral", "Lateral cortex · front at left", "dk", (-5, 625, 780, 1180)), ("medial", "Medial cortex · front at right", "dk", (800, 655, 1610, 1175)),
         ("coronal", "Coronal slice · striatum, pallidum, amygdala", "aseg", (-5, -5, 1345, 1240)),
         ("sagittal", "Midsagittal slice · brainstem, cerebellum · front at right", "aseg", (1525, 60, 2900, 1240))]


def _subpaths(d):
    import re
    out = []
    for sp in re.split(r"(?=M )", d.strip()):
        nums = [float(x) for x in re.findall(r"-?\d+\.?\d*", sp)]
        if len(nums) >= 2:
            xs, ys = nums[0::2], nums[1::2]
            out.append((sp.strip(), (min(xs), max(xs), min(ys), max(ys)), (sum(xs) / len(xs), sum(ys) / len(ys))))
    return out


def build_anat():
    """Four atlas views with background parcels, coloured sites, and non-overlapping margin labels."""
    site_of = {}
    for site, parts in ANAT_MAP.items():
        for a, f, *c in parts:
            site_of.setdefault((a, f), []).append((site, c[0] if c else None))
    views = []
    for key, title, atlas, (x0, y0, x1, y1) in VIEWS:
        inside = lambda b: b[0] >= x0 - 5 and b[1] <= x1 + 5 and b[2] >= y0 - 5 and b[3] <= y1 + 5  # noqa: E731
        base, parcels, sites = [], [], {}
        for f in sorted((ATLAS / atlas).iterdir()):
            if f.name == "README.md":
                continue
            for d, bb, cen in _subpaths(f.read_text()):
                if not inside(bb):
                    continue
                if f.name in ("lateral_left", "medial_left", "Coronal", "Sagittal"):
                    base.append(d)
                    continue
                hits = site_of.get((atlas, f.name))
                if hits:
                    if any(clip for _, clip in hits):  # whole structure as neutral background, coloured bands on top
                        parcels.append({"site": None, "d": d, "name": f.name})
                    for site, clip in hits:
                        item = {"site": site, "d": d}
                        cy = cen[1]
                        if clip:
                            lo, hi = clip
                            ylo, yhi = (lo if lo is not None else bb[2] - 10), (hi if hi is not None else bb[3] + 10)
                            item["clip"] = [bb[0] - 10, ylo, bb[1] + 10, yhi]
                            cy = (max(ylo, bb[2]) + min(yhi, bb[3])) / 2
                        parcels.append(item)
                        sites.setdefault(site, []).append((cen[0], cy))
                else:
                    parcels.append({"site": None, "d": d, "name": f.name})
        w, h = x1 - x0, y1 - y0
        fs = max(w, 900) / 38
        labels = []
        for side in ("left", "right"):
            items = []
            for site, pts in sites.items():
                cx = sum(p[0] for p in pts) / len(pts); cy = sum(p[1] for p in pts) / len(pts)
                # one anchor per site: the left-hemisphere copy for bilateral coronal structures
                if key == "coronal":
                    left = min(pts, key=lambda p: p[0]); right = max(pts, key=lambda p: p[0]); cx, cy = (left if side == "left" else right)
                    want = "left" if site in ("putamen", "amygdala") else "right"
                else:
                    want = "left" if cx < x0 + w / 2 else "right"
                if want == side:
                    items.append([site, cx, cy])
            items.sort(key=lambda t: t[2])
            ly = y0 + fs
            for it in items:
                ly = max(ly, it[2]); labels.append({"site": it[0], "ax": it[1], "ay": it[2], "lx": x0 - fs * .6 if side == "left" else x1 + fs * .6,
                                                   "ly": ly, "anchor": "end" if side == "left" else "start", "text": ANAT_LABEL.get(it[0], it[0])})
                ly += fs * 2.6
        braak = {"prefrontal cortex": 5, "motor cortex": 6, "cingulate cortex": 4, "visual cortex": 6, "temporal cortex": 4, "amygdala": 4,
                 "substantia nigra": 3, "vagal nucleus (DMNX)": 1}
        need = lambda side: max([len(L["text"]) + (5 if L["site"] in braak else 0) for L in labels if L["anchor"] == side] or [0]) * fs * .56 + fs  # noqa: E731
        ml, mr = need("end"), need("start")
        views.append({"key": key, "title": title, "vb": [x0 - ml, y0 - fs, w + ml + mr, h + 3 * fs], "fs": round(fs, 1),
                      "base": base, "parcels": parcels, "labels": labels})
    return views


def recs(df, digits=4):
    df = df.copy()
    for c in df.select_dtypes("number"):
        df[c] = df[c].round(digits)
    return df.replace({np.nan: None}).to_dict("records")


OPTIONAL = {"consensus_genes_PDany.tsv": ["gene", "cell_class", "est", "se", "fdr", "k", "n_sig_same_dir", "frac_dir"],
            "gsea_PDany_vs_control.tsv.gz": ["cell_class", "library", "term", "nes", "fdr", "n_genes", "scope"]}


def tsv(name, **kw):
    if not (BIO / name).exists() and BIO != RES / "biology" and name in OPTIONAL:  # public tables still being computed
        print(f"WARNING: {BIO / name} not built yet; its section stays empty", file=sys.stderr)
        return pd.DataFrame(columns=OPTIONAL[name])
    return pd.read_csv(BIO / name, sep="\t", low_memory=False, **kw)


def figure_b64(name, width=1400):
    from PIL import Image
    im = Image.open(ROOT / "figures" / name).convert("RGB")
    if im.width > width:
        im = im.resize((width, int(im.height * width / im.width)))
    buf = BytesIO()
    im.save(buf, "JPEG", quality=82, optimize=True)
    return "data:image/jpeg;base64," + base64.b64encode(buf.getvalue()).decode()


def store():
    con = sqlite3.connect(DB)
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
    global BG_RECORDS
    BG_RECORDS = bgt.replace({np.nan: None}).to_dict("records")
    sig = q("""SELECT cs.claim_id c, sg.gene g FROM claim_signature cs
               JOIN signature_gene sg ON cs.signature_id = sg.signature_id""").dropna().drop_duplicates()
    papers = q("SELECT paper_id id, first_author, year, title, doi, journal_or_preprint journal FROM source "
               "WHERE paper_id IN (SELECT DISTINCT paper FROM claim)")
    info = q("SELECT * FROM build_info").iloc[0].to_dict()
    con.close()
    named = {}
    for _, r in claims.iterrows():
        for g in filter(None, (r.genes or "").split(";")):
            named.setdefault(g, set()).add(r.id)
    in_sig = sig.groupby("g").c.apply(set).to_dict()
    gene_claims = {g: [sorted(named.get(g, set())), sorted(in_sig.get(g, set()) - named.get(g, set()))]
                   for g in set(named) | set(in_sig)}
    return claims, ds, ev, rel, cs, papers, info, gene_claims


def question_form():
    """GitHub issue form for questions: of the public repository (config [release] public_repository), else of the
    repository in CITATION.cff."""
    import re
    if cfg.PUBLIC_REPOSITORY:
        return cfg.PUBLIC_REPOSITORY.rstrip("/") + "/issues/new?template=question.yml"
    m = re.search(r'^repository-code:\s*"?([^"\n]+)"?', (ROOT / "CITATION.cff").read_text(), re.M)
    return m.group(1).rstrip("/") + "/issues/new?template=question.yml" if m else None


def main():
    global DB, BIO, PSY
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--public", action="store_true", help="public atlas (see the module docstring)")
    ap.add_argument("--db", help="claim store (default: the configured store; required with --public)")
    a = ap.parse_args()
    if a.db:
        DB = Path(a.db)
    if a.public:
        if not a.db:
            sys.exit("--public needs --db: the claim_store.sqlite of a public release (scripts/59_release.py)")
        BIO, PSY = RES / "public" / "biology", RES / "public" / "psychad"
        for k in cfg.PRIVATE_MODULES:
            PROGRAMMES.pop(k, None)
    claims, ds, ev, rel, cs, papers, info, gene_claims = store()
    src = pd.read_sql("SELECT paper_id, first_author, year FROM source", sqlite3.connect(DB)).set_index("paper_id")
    def label(d):
        a = f"{src.first_author.get(d.dataset, '')} {int(float(src.year.get(d.dataset)))}" if d.dataset in src.index and src.year.get(d.dataset) else d.dataset
        return a if d.cohort == d.dataset else f"{a} · {d.cohort}"
    ds["label"] = ds.apply(label, axis=1)

    comp = tsv("composition_meta.tsv")[["contrast", "cell_class", "denominator", "scope", "est", "se", "p", "fdr",
                                        "k", "I2", "datasets"]]
    prog = tsv("programme_meta.tsv")[["programme", "cell_class", "test", "est", "se", "p", "fdr", "k", "I2",
                                      "datasets"]]
    region = tsv("region_map.tsv")[["programme", "cell_class", "site", "est", "se", "p", "k", "donors", "datasets"]]
    spec = pd.read_csv(PSY / "specificity.tsv", sep="\t")[
        ["programme", "cell_class", "group", "effect_sd", "se_sd", "p", "fdr", "n_group", "n_ctrl",
         "share_above_ctrl95"]]
    gen = tsv("genetics_meta.tsv")[["set_type", "cell_class", "k", "z", "p", "share_p05", "datasets_p05"]]
    het = tsv("donor_heterogeneity.tsv")[["programme", "cell_class", "share_above", "share_below", "n_cases",
                                          "log_sd_ratio", "p"]]
    # the old theme summary (evidence_map) counted private claims and tests; the page now groups themes from the claims
    themes = pd.DataFrame() if a.public else tsv("evidence_map.tsv")
    genes = tsv("consensus_genes_PDany.tsv")[["gene", "cell_class", "est", "se", "fdr", "k", "n_sig_same_dir",
                                              "frac_dir"]]
    gsea = tsv("gsea_PDany_vs_control.tsv.gz")
    gsea = gsea[(gsea.scope == "all regions") & gsea.library.isin(["Hallmark", "Reactome"]) & (gsea.fdr < 0.05)]
    gsea = (pd.concat([gsea[gsea.nes > 0].sort_values("nes", ascending=False).groupby("cell_class").head(10),
                       gsea[gsea.nes < 0].sort_values("nes").groupby("cell_class").head(10)])
            [["cell_class", "library", "term", "nes", "fdr", "n_genes"]])
    if a.public:  # safeguard: no row of a private paper's gene module, whatever the input tables hold
        prog, region, spec, het = (t[~t.programme.isin(cfg.PRIVATE_MODULES)] for t in (prog, region, spec, het))
    gsea["term"] = (gsea.term.str.replace(r"^(HALLMARK|REACTOME)_", "", regex=True).str.replace("_", " ")
                    .str.lower().str.capitalize())

    # headline numbers
    tested = claims.verdict.fillna("Untested").ne("Untested")
    held = claims.verdict.fillna("").str.startswith("Supported")
    strong = held & claims.certainty.isin(["High", "Moderate"])
    regions = sorted({r for x in ds.regions for r in x.split(";")} - {"other"})
    findings = claims.origin.eq("finding")
    stats = dict(findings=int(findings.sum()), papers=int(claims.paper[~findings].nunique()), datasets=len(ds), donors=int(ds.n_donors.sum()),
                 regions=len(regions), claims=int((~findings).sum()), tested=int((tested & ~findings).sum()), held=int((held & ~findings).sum()), strong=int((strong & ~findings).sum()),
                 nuclei=None if a.public else "10.9 million", nuclei_source=None if a.public else "manuscript v0.8 abstract")

    data = dict(
        claims=recs(claims), datasets=recs(ds), evidence=recs(ev, 3), relations=recs(rel), statements=recs(cs),
        background=BG_RECORDS,
        papers=recs(papers), gene_claims=gene_claims, info={k: str(v) for k, v in info.items()}, stats=stats,
        comp=recs(comp), prog=recs(prog), region=recs(region), spec=recs(spec), gen=recs(gen), het=recs(het),
        themes=recs(themes), genes=recs(genes, 3), gsea=recs(gsea, 3),
        programmes={k: dict(group=g, name=n) for k, (g, n) in PROGRAMMES.items()}, cell_order=CELL_ORDER,
        sites={k: dict(x=x, y=y, braak=b, regions=SITE_REGIONS[k], dx=SITE_LAB[k][0], dy=SITE_LAB[k][1], anchor=SITE_LAB[k][2])
               for k, (x, y, b) in SITES.items()},
        model_figure=None if a.public else figure_b64("B9_model.png"), anat=build_anat(),
        public=bool(a.public), links={"question_form": question_form()} if a.public else {},
        kg=kg.build(DB).compact(kg.similarity_links(DB, cache=cfg.STORE_DB.parent / "claim_embeddings.npz"),
                                kg.topic_clusters(DB, map_dir=RES / "public" / "claim_map" if a.public else None)))
    tpl = (ROOT / "scripts" / "atlas_template.html").read_text()
    out = ROOT / "reports" / ("pd_atlas_public.html" if a.public else "pd_atlas.html")
    html = tpl.replace("/*__DATA__*/null", json.dumps(data, separators=(",", ":"), default=str))
    out.write_text(html)
    if a.public:
        print(f"{out} ({out.stat().st_size / 1e6:.1f} MB): {len(claims)} claims, {len(ds)} datasets, {len(ev)} tests")
        return
    # claude.ai Artifact version: the publisher adds the doctype/head/body skeleton, so ship title + style + body only
    head, body = html.split("<body>", 1)
    art = head[head.index("<title>"):head.index("</head>")] + body.rsplit("</body>", 1)[0]
    (ROOT / "reports" / "pd_atlas.artifact.html").write_text(art)
    print(f"{out} ({out.stat().st_size / 1e6:.1f} MB): {len(claims)} claims, {len(ds)} datasets, {len(ev)} tests, "
          f"{len(prog)} programme effects, {len(genes)} gene effects, {len(gsea)} pathways")


if __name__ == "__main__":
    main()

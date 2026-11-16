"""Disease-programme map: preranked GSEA on the consensus (meta-analysis) gene effects per cell type and region.

Usage: python -I scripts/24_pathways.py <contrast> [cell class]   (contrast: PDany_vs_control or braak)
With a cell class, only that class is run and written to gsea_parts/<contrast>__<class>.tsv.gz, so classes can run as
parallel SLURM jobs; without one, all parts present are merged into gsea_<contrast>.tsv.gz.
Ranks: meta z of results/biology/de_meta_<contrast>.tsv.gz per cell class x scope (all regions, midbrain / SN,
striatum, cortex). Gene sets: MSigDB 2024.1 Hallmark and Reactome for every scope, GO biological process for the
all-region rankings only (15-500 genes).
GSEA: gseapy.prerank, 1000 permutations, seed 0. Writes results/biology/gsea_<contrast>.tsv.gz.
"""
from pathlib import Path
import sys
import pandas as pd
import gseapy as gp

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

BASE = cfg.WORK
GS = BASE / "raw" / "genesets"
# --public: rank the public gene meta-analysis (23 --public) and write next to it, in results/public/biology
PUBLIC = "--public" in sys.argv
OUT = cfg.RESULTS / "public" / "biology" if PUBLIC else BASE / "results" / "biology"
LIBS = {"Hallmark": "h.all.v2024.1.Hs.symbols.gmt", "Reactome": "c2.cp.reactome.v2024.1.Hs.symbols.gmt",
        "GO_BP": "c5.go.bp.v2024.1.Hs.symbols.gmt"}
MIN_GENES = 3000


def merge(contrast):
    parts = sorted((OUT / "gsea_parts").glob(f"{contrast}__*.tsv.gz"))
    out = pd.concat([pd.read_csv(f, sep="\t") for f in parts], ignore_index=True)
    out.to_csv(OUT / f"gsea_{contrast}.tsv.gz", sep="\t", index=False)
    h = out[(out.library == "Hallmark") & (out.scope == "all regions") & (out.fdr < 0.05)]
    print(len(parts), "parts;", h.pivot_table(index="term", columns="cell_class", values="nes").round(1).to_string())


def main(contrast, only=None):
    m = pd.read_csv(OUT / f"de_meta_{contrast}.tsv.gz", sep="\t")
    if only:
        m = m[m.cell_class == only]
    res = []
    for (cls, scope), d in m.groupby(["cell_class", "scope"]):
        rnk = d.dropna(subset=["z"]).drop_duplicates("gene").set_index("gene").z.sort_values(ascending=False)
        if len(rnk) < MIN_GENES:
            continue
        for lib, f in LIBS.items():
            if lib == "GO_BP" and scope != "all regions":  # GO BP only for all-region rankings (runtime)
                continue
            r = gp.prerank(rnk=rnk, gene_sets=str(GS / f), min_size=15, max_size=500, permutation_num=1000,
                           seed=0, threads=4, outdir=None, verbose=False).res2d
            r = r.rename(columns={"Term": "term", "NES": "nes", "NOM p-val": "p", "FDR q-val": "fdr",
                                  "Lead_genes": "leading_edge", "Tag %": "tag"})
            res.append(r[["term", "nes", "p", "fdr", "leading_edge"]].assign(cell_class=cls, scope=scope, library=lib,
                                                                         n_genes=len(rnk)))
        print(cls, scope, flush=True)
    if not res:
        return
    out = pd.concat(res, ignore_index=True)
    for c in ["nes", "p", "fdr"]:
        out[c] = pd.to_numeric(out[c], errors="coerce")
    (OUT / "gsea_parts").mkdir(exist_ok=True)
    out.to_csv(OUT / "gsea_parts" / f"{contrast}__{only.replace(' ', '_')}.tsv.gz", sep="\t", index=False)


if __name__ == "__main__":
    args = [a for a in sys.argv[1:] if a != "--public"]
    con = args[0] if args else "PDany_vs_control"
    if len(args) > 1:
        main(con, args[1])
    else:
        merge(con)

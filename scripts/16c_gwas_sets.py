"""Per-dataset gene sets for PD GWAS enrichment, and the MAGMA competitive gene-set test on them.

Usage: python -I scripts/16c_gwas_sets.py <dataset>
Sets (gene symbols mapped to Entrez via the MAGMA NCBI37.3 gene locations):
  spec::<class>               top 10% most cell-class-specific genes in control units (share of the gene's mean CPM
                              across classes; genes with mean CPM >= 1 in the class), cf. Bryois 2020 / Kamath 2022
  deup::<contrast>::<class>   genes up in cases (edgeR FDR < 0.05; if fewer than 50, the top 200 by p with logFC > 0)
Writes results/gwas/sets/<dataset>.sets (MAGMA set-annot), <dataset>.gsa.out (MAGMA), <dataset>.spec.tsv.
Also writes the familial-PD-gene overlap (hypergeometric) per spec set to <dataset>.familial.tsv.
"""
from pathlib import Path
import subprocess
import sys
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.harmonise import OUT as H  # noqa: E402
from pdrep.contrasts import unit_groups  # noqa: E402

from pdrep import config as cfg  # noqa: E402

BASE = cfg.WORK
G = BASE / "raw" / "gwas"
GW = BASE / "results" / "gwas"
DE = BASE / "results" / "de_cache"
MAGMA = G / "magma_bin" / "x" / "magma"
# Mendelian / familial PD genes (MDS gene list and Kamath 2022 familial set, current HGNC symbols)
FAMILIAL = ["SNCA", "LRRK2", "PRKN", "PINK1", "PARK7", "VPS35", "GBA1", "ATP13A2", "PLA2G6", "FBXO7", "DNAJC6",
            "SYNJ1", "VPS13C", "CHCHD2", "DCTN1", "GCH1", "TH", "SPR", "UCHL1", "GIGYF2", "HTRA2", "EIF4G1", "DNAJC13",
            "TMEM230", "RAB39B", "LRP10", "POLG", "ATXN2", "ATXN3", "MAPT"]
CONTRASTS = ["PDany_vs_control", "PDD_vs_control", "LBDany_vs_control", "LRRK2_vs_control", "GBA_vs_control"]


def entrez_map():
    loc = pd.read_csv(G / "magma_ref" / "x" / "NCBI37.3.gene.loc", sep="\t", header=None,
                      names=["entrez", "chr", "start", "end", "strand", "symbol"])
    return dict(zip(loc.symbol.str.upper(), loc.entrez.astype(str)))


def main(ds):
    out = GW / "sets"
    out.mkdir(parents=True, exist_ok=True)
    e = entrez_map()
    pb = ad.read_h5ad(H / "pseudobulk" / f"{ds}.h5ad")
    if "enrichment" in pb.obs:
        pb = pb[~pb.obs.enrichment.astype(str).str.lower().eq("sorted")].copy()
    g = unit_groups(pb.obs.assign(disease=pb.obs.disease.astype(str)))
    ctrl = pb[(pb.obs.disease.astype(str) == "control").values]
    if ctrl.n_obs == 0:
        ctrl = pb
    X = ctrl.X.toarray() if hasattr(ctrl.X, "toarray") else ctrl.X
    cpm = X / X.sum(1, keepdims=True) * 1e6
    m = pd.DataFrame(cpm, columns=ctrl.var_names).groupby(ctrl.obs.cell_class.astype(str).values).mean()
    m = m[[c for c in m.columns]]
    spec = m / m.sum(0).replace(0, np.nan)
    sets, fam = {}, []
    universe = [x for x in m.columns if x.upper() in e]
    famset = {x.upper() for x in FAMILIAL}
    for cls in spec.index:
        ok = m.loc[cls] >= 1
        s = spec.loc[cls][ok].dropna().sort_values(ascending=False)
        top = list(s.index[: max(50, int(0.10 * ok.sum()))])
        sets[f"spec::{cls}"] = top
        k = sum(x.upper() in famset for x in top)
        K = sum(x.upper() in famset for x in universe)
        p = stats.hypergeom.sf(k - 1, len(universe), K, len(top)) if k else 1.0
        fam.append({"set": f"spec::{cls}", "cell_class": cls, "n_genes": len(top), "familial_in_set": k,
                    "familial_in_universe": K, "p": p,
                    "genes": ";".join(sorted(x for x in top if x.upper() in famset))})
    for f in (DE / ds / "de").glob("*__*.tsv.gz"):
        reg_cls, con = f.name[: -len(".tsv.gz")].rsplit("__", 1)
        if con not in CONTRASTS:
            continue
        cls = reg_cls.split("__", 1)[1].replace("_", " ")
        t = pd.read_csv(f, sep="\t")
        up = t[(t.padj < 0.05) & (t.logFC > 0)].gene.tolist()
        if len(up) < 50:
            up = t[t.logFC > 0].nsmallest(200, "p").gene.tolist()
        sets.setdefault(f"deup::{con}::{cls}", [])
        sets[f"deup::{con}::{cls}"] = list(dict.fromkeys(sets[f"deup::{con}::{cls}"] + up))  # union over regions
    with open(out / f"{ds}.sets", "w") as fh:
        for name, genes in sets.items():
            ids = [e[x.upper()] for x in genes if x.upper() in e]
            if len(ids) >= 10:
                fh.write(name.replace(" ", "_") + " " + " ".join(dict.fromkeys(ids)) + "\n")
    pd.DataFrame(fam).to_csv(out / f"{ds}.familial.tsv", sep="\t", index=False)
    spec.T.to_csv(out / f"{ds}.spec.tsv.gz", sep="\t")
    subprocess.run([str(MAGMA), "--gene-results", str(GW / "nalls2019.genes.raw"), "--set-annot", str(out / f"{ds}.sets"),
                    "--out", str(out / ds)], check=True)
    print(ds, len(sets), "sets")


if __name__ == "__main__":
    main(sys.argv[1])

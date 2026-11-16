"""Disease programmes per donor: what do they track (Lewy pathology, clinical PD, dementia, AD, or post-mortem factors)?

Usage: python -I scripts/27_programmes.py
For every dataset x region x cell class pseudobulk (results/de_cache/<ds>/input), genes are log2(CPM+1) and z-scored
across donor units; a programme score is the mean z of its genes (>= 5 genes present; 2 for the cold-inducible pair).
Tests per dataset x region x class: Hedges' g for each group contrast (pdrep.contrasts) and Spearman rho (as Fisher z)
with Lewy-body stage (controls = 0), post-mortem interval and age. Regions of a dataset are combined by fixed effect,
datasets by DerSimonian-Laird. Writes results/biology/programme_effects.tsv and programme_meta.tsv.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.contrasts import braak_value  # noqa: E402
from pdrep.meta import fixed, dersimonian_laird, bh  # noqa: E402

from pdrep import config as cfg  # noqa: E402

BASE = cfg.WORK
DE, GS, OUT = BASE / "results" / "de_cache", BASE / "raw" / "genesets", BASE / "results" / "biology"
CONTRASTS = ["PDany_vs_control", "PD_vs_control", "PDD_vs_control", "DLB_vs_control", "iLBD_vs_control",
             "LBearly_vs_control", "LBmid_vs_control", "LBlate_vs_control", "AD_vs_control"]
MSIGDB = {"heat response": ("GOBP_CELLULAR_RESPONSE_TO_HEAT", "c5.go.bp"),
          "chaperone folding": ("GOBP_CHAPERONE_MEDIATED_PROTEIN_FOLDING", "c5.go.bp"),
          "unfolded protein response": ("HALLMARK_UNFOLDED_PROTEIN_RESPONSE", "h.all"),
          "hypoxia": ("HALLMARK_HYPOXIA", "h.all"),
          "inflammatory response": ("HALLMARK_INFLAMMATORY_RESPONSE", "h.all"),
          "TNFa / NF-kB": ("HALLMARK_TNFA_SIGNALING_VIA_NFKB", "h.all"),
          "interferon gamma": ("HALLMARK_INTERFERON_GAMMA_RESPONSE", "h.all"),
          "oxidative phosphorylation": ("HALLMARK_OXIDATIVE_PHOSPHORYLATION", "h.all")}
STATES = {  # published state signatures (claims/signatures_all.tsv.gz) and the classes they are scored in
    "UPR astrocytes (Barba-Reyes)": ("barbareyes2026__astro_UPR_module1", ["Astrocytes"]),
    "reactive astrocytes (Barba-Reyes)": ("barbareyes2026__astro_Reactive_module3__preprint", ["Astrocytes"]),
    "UPR microglia (Barba-Reyes)": ("barbareyes2026__micro_UPR_module9", ["Microglia"]),
    "reactive microglia (Barba-Reyes)": ("barbareyes2026__micro_Reactive_module5", ["Microglia"]),
    "DAM microglia (Barba-Reyes)": ("barbareyes2026__DAM_markers", ["Microglia"]),
    "GPNMB-high microglia (Smajić)": ("smajic2022__microglia_GPNMBhigh_branch", ["Microglia"]),
    "AR2 astrocytes (Barba-Reyes)": ("barbareyes2026__AR2_markers", ["Astrocytes"]),
    "heat-shock oligodendroglia (Barba-Reyes)": ("barbareyes2025__HSPR_module_xenium", ["Oligodendrocytes"]),
    "SOX6/AGTR1 DA subtype (Kamath)": ("kamath2022__DA_SOX6_AGTR1_markers", ["DA neurons"]),
    "HSP90AA1-high microglia (Smajić)": ("smajic2022__microglia_HSP90AA1high_branch", ["Microglia"]),
    "C3/CD44 astrocytes, Astrocytes1 (Martirosyan)": ("martirosyan2024__astrocytes1_markers", ["Astrocytes"]),
    "UPR astrocytes, Astrocytes3 (Martirosyan)": ("martirosyan2024__astrocytes3_markers", ["Astrocytes"]),
    "GFAP/MT astrocytes, Astrocytes4 (Martirosyan)": ("martirosyan2024__astrocytes4_markers", ["Astrocytes"]),
    "stress oligodendrocytes, Oligos5 (Martirosyan)": ("martirosyan2024__oligos5_markers", ["Oligodendrocytes"]),
    "activated astrocytes, Ast-1 (Xu)": ("xu2023__Ast1_markers", ["Astrocytes"]),
    "activated microglia, Micr-1 (Xu)": ("xu2023__Micr1_markers", ["Microglia"]),
    "activated myeloid cells (Ma)": ("ma2025__SN_Myeloid_Activated_markers", ["Microglia"]),
    "disease-associated oligodendrocytes (Adams)": ("adams2024__ODC_disease_assoc_markers", ["Oligodendrocytes"]),
}
CLASSES = ["Astrocytes", "Microglia", "Oligodendrocytes", "OPCs", "Endothelial", "Excitatory neurons",
           "Inhibitory neurons", "DA neurons", "MSNs", "ALL"]


def gmt(name, lib):
    for line in open(GS / f"{lib}.v2024.1.Hs.symbols.gmt"):
        f = line.rstrip("\n").split("\t")
        if f[0] == name:
            return [g.upper() for g in f[2:]]
    raise KeyError(name)


def programmes():
    p = {k: (gmt(n, lib), CLASSES) for k, (n, lib) in MSIGDB.items()}
    p["cold-inducible (CIRBP, RBM3)"] = (["CIRBP", "RBM3"], CLASSES)
    sig = pd.read_csv(ROOT / "claims" / "signatures_all.tsv.gz", sep="\t", dtype=str, usecols=["signature_id", "gene"])
    for k, (sid, cls) in STATES.items():
        p[k] = (sig[sig.signature_id == sid].gene.str.upper().unique().tolist(), cls)
    return p


def hedges(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    n1, n2 = len(a), len(b)
    if n1 < 3 or n2 < 3:
        return np.nan, np.nan
    sp = np.sqrt(((n1 - 1) * a.var(ddof=1) + (n2 - 1) * b.var(ddof=1)) / (n1 + n2 - 2))
    if sp == 0:
        return np.nan, np.nan
    j = 1 - 3 / (4 * (n1 + n2) - 9)
    g = j * (a.mean() - b.mean()) / sp
    return g, np.sqrt((n1 + n2) / (n1 * n2) + g ** 2 / (2 * (n1 + n2)))


def fisher(x, y):
    ok = np.isfinite(x) & np.isfinite(y)
    if ok.sum() < 6 or np.unique(y[ok]).size < 3:
        return np.nan, np.nan
    rho = stats.spearmanr(x[ok], y[ok])[0]
    return np.arctanh(np.clip(rho, -0.999, 0.999)), 1 / np.sqrt(ok.sum() - 3)


def pool(eff):
    """Per dataset x region effects -> per programme x cell class x test: one estimate per dataset (regions by fixed
    effect), DerSimonian-Laird across datasets; BH FDR per test across programmes and cell classes."""
    res = []
    for (prog, cls, test), d in eff.groupby(["programme", "cell_class", "test"]):
        per = d.groupby("dataset").apply(lambda x: pd.Series(fixed(x.est, x.se), index=["est", "se"]))
        r = dersimonian_laird(per.est, per.se)
        res.append({"programme": prog, "cell_class": cls, "test": test, **r, "datasets": ";".join(per.dropna().index)})
    meta = pd.DataFrame(res)
    meta["fdr"] = meta.groupby("test").p.transform(lambda p: pd.Series(bh(p.values), index=p.index))
    return meta


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    progs = programmes()
    rows = []
    for ds_dir in sorted(DE.iterdir()):
        units_f = ds_dir / "input" / "units.tsv"
        if not units_f.exists():
            continue
        units = pd.read_csv(units_f, sep="\t", dtype=str).set_index("unit")
        for f in sorted((ds_dir / "input").glob("*__*.csv.gz")):
            region, cls_f = f.name[: -len(".csv.gz")].split("__", 1)
            cls = cls_f.replace("_", " ").replace("Neurons other", "Neurons (other)")
            if cls not in CLASSES:
                continue
            x = pd.read_csv(f, index_col=0)
            x.index = x.index.str.upper()
            lc = np.log2(x.div(x.sum(0), axis=1) * 1e6 + 1)
            z = lc.sub(lc.mean(1), axis=0).div(lc.std(1).replace(0, np.nan), axis=0)
            u = units.loc[z.columns]
            stage = braak_value(u.assign(braak_lb=pd.to_numeric(u.braak_lb, errors="coerce"))).values
            for pname, (genes, classes) in progs.items():
                if cls not in classes:
                    continue
                g = [x_ for x_ in genes if x_ in z.index]
                if len(g) < (2 if len(genes) == 2 else 5):
                    continue
                score = z.loc[g].mean(0).values
                base = {"dataset": ds_dir.name, "region": region, "cell_class": cls, "programme": pname, "n_genes": len(g)}
                for con in CONTRASTS:
                    if con in u:
                        grp = u[con].fillna("").values
                        est, se = hedges(score[grp == "case"], score[grp == "ref"])
                        if np.isfinite(est):
                            rows.append({**base, "test": con, "est": est, "se": se, "n": int((grp != "").sum())})
                for test, cov in [("Lewy stage", stage), ("post-mortem interval", pd.to_numeric(u.pmi, errors="coerce").values),
                                  ("age", pd.to_numeric(u.age, errors="coerce").values)]:
                    est, se = fisher(score, cov)
                    if np.isfinite(est):
                        rows.append({**base, "test": test, "est": est, "se": se, "n": int(np.isfinite(cov).sum())})
        print(ds_dir.name, flush=True)
    eff = pd.DataFrame(rows)
    eff.to_csv(OUT / "programme_effects.tsv", sep="\t", index=False)
    meta = pool(eff)
    meta.to_csv(OUT / "programme_meta.tsv", sep="\t", index=False)
    show = meta[(meta.cell_class.isin(["ALL", "Astrocytes", "Microglia", "Oligodendrocytes"])) & (meta.k >= 2)]
    print(show.pivot_table(index=["programme", "cell_class"], columns="test", values="z").round(1).to_string())


if __name__ == "__main__":
    main()

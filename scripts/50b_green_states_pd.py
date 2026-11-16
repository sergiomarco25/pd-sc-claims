"""Do Lewy-body brains move along the AD trajectory of Green et al. 2024, the alternative ageing trajectory, or neither?

Usage: python -I scripts/50b_green_states_pd.py [pd|pd_nocore|psychad]
Signatures: for each state of Green et al. 2024 that defines a trajectory (prAD: Mic.12, Mic.13, Ast.10, Oli.7, Mic.11,
Ast.9; ABA: Ast.5, Mic.10, OPC.3, Oli.8, Mic.6; homeostatic C1: Ast.1, Oli.4 and others), the top 25 one-vs-rest
markers of Supp. Table 2 (adjusted p < 0.05, expressed in >= 20% of the state's nuclei, log2FC < 10 to drop the rare
artefactual fold changes) are scored in the matching cell class.
  pd:      every dataset of the meta-analysis, as scripts/27_programmes.py (Hedges' g per contrast, DerSimonian-Laird
           across datasets); writes results/biology/green2024_pd/programme_effects.tsv and programme_meta.tsv.
  pd_nocore: as pd, with heat-shock/chaperone genes and the shared PD core genes removed from every signature
           (results/biology/green2024_pd_nocore).
  psychad: all PsychAD diagnoses versus controls, as scripts/41b_psychad_specificity.py; plus, across donors, the
           Spearman correlation of each state score with the donor stress index and the cold-inducible score.
           Writes results/psychad/green2024_specificity.tsv and green2024_stress_correlation.tsv.
"""
from pathlib import Path
import importlib.util
import sys
import warnings
import re
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parents[1]
warnings.simplefilter("ignore")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

RAW = cfg.RAW / "green2024"
RES = cfg.RESULTS
CLS = {"Mic": "Microglia", "Ast": "Astrocytes", "Oli": "Oligodendrocytes", "OPC": "OPCs", "End": "Endothelial",
       "Exc": "Excitatory neurons", "Inh": "Inhibitory neurons"}
STATES = {"Mic.12": "prAD early", "Mic.13": "prAD early", "Mic.11": "prAD late", "Ast.9": "prAD late",
          "Ast.10": "prAD late", "Oli.7": "prAD late", "End.5": "prAD late", "OPC.1": "prAD mid",
          "Ast.5": "ABA", "Mic.10": "ABA", "Mic.6": "ABA", "Mic.15": "ABA", "OPC.3": "ABA", "Oli.8": "ABA",
          "Exc.1": "ABA", "Ast.1": "homeostatic", "Mic.2": "homeostatic", "Oli.4": "ABA", "Mic.14": "homeostatic",
          "Ast.7": "homeostatic", "Inh.6": "homeostatic", "Inh.16": "prAD neuronal"}


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


CORE = re.compile(r"^(HSPA|HSPH|HSP90|HSPB|HSPD|HSPE|DNAJ|BAG3$|CRYAB$|CHORDC1$|STIP1$|AHSA1$|PTGES3$|FKBP4$|CACYBP$|ST13$|"
                  r"SERPINH1$|SLC38A2$|NEAT1$|MALAT1$|CIRBP$|RBM3$)")


def green_programmes(n=25, drop_core=False):
    """Top-n markers per state; drop_core removes heat-shock/chaperone genes and the shared PD core (SLC38A2, NEAT1,
    MALAT1, CIRBP, RBM3) before ranking, so the score asks whether the rest of the state moves too."""
    d = pd.read_csv(RAW / "green2024_suppT2_DEGs_one_vs_rest.csv.gz")
    if drop_core:
        d = d[~d.gene.astype(str).str.upper().str.match(CORE)]
    d = d[(d.p_val_adj < 0.05) & (d.pct_1 if "pct_1" in d else d["pct.1"]).ge(0.2) & d.avg_log2FC.between(0.25, 10)]
    out = {}
    for s, traj in STATES.items():
        g = d[d.state == s].sort_values("avg_log2FC", ascending=False).gene.astype(str).str.upper().drop_duplicates()
        out[f"{s} ({traj})"] = (g.head(n).tolist(), [CLS[s.split(".")[0]]])
    return out


def run_pd(drop_core=False):
    prog = load("prog", ROOT / "scripts" / "27_programmes.py")
    gp = green_programmes(drop_core=drop_core)
    prog.programmes = lambda: gp
    prog.OUT = RES / "biology" / ("green2024_pd_nocore" if drop_core else "green2024_pd")
    prog.main()


def run_psychad():
    sp = load("spec41", ROOT / "scripts" / "41b_psychad_specificity.py")
    base = sp.prog.programmes()
    keep = {k: base[k] for k in ["heat response", "chaperone folding", "cold-inducible (CIRBP, RBM3)"]}
    sp.prog.programmes = lambda: {**keep, **green_programmes()}
    meta = pd.read_csv(sp.P / "donors.tsv", sep="\t")
    s = sp.scores().merge(meta[["donor", "group", "age", "sex"]], on="donor")
    gl = s[s.programme.isin(["heat response", "chaperone folding"])
           & s.cell_class.isin(["Astrocytes", "Microglia", "Oligodendrocytes", "OPCs"])].copy()
    ref = gl[gl.group == "control"].groupby(["programme", "cell_class"]).score.agg(["mean", "std"])
    gl = gl.join(ref, on=["programme", "cell_class"])
    gl["zs"] = (gl.score - gl["mean"]) / gl["std"]
    idx = gl.groupby("donor").zs.mean().rename("stress_index")
    cold = s[(s.programme == "cold-inducible (CIRBP, RBM3)") & s.cell_class.isin(["Astrocytes", "Oligodendrocytes"])]
    cold = cold.groupby("donor").score.mean().rename("cold_inducible")
    rows, cor = [], []
    for (p, c), d in s.groupby(["programme", "cell_class"]):
        if not p[0].isupper() or "." not in p:
            continue
        for g in sp.GROUPS:
            r = sp.fit(d, g)
            if r:
                rows.append({"programme": p, "cell_class": c, "group": g, **r})
        x = d.set_index("donor").join(idx).join(cold)
        for grp_name, sub in [("all donors", x), ("controls", x[x.group == "control"]),
                              ("AD", x[x.group == "Alzheimer's disease"]),
                              ("Lewy-body (with or without AD)", x[x.group.str.startswith("Lewy")])]:
            for col in ["stress_index", "cold_inducible"]:
                ok = sub[["score", col]].dropna()
                if len(ok) >= 20:
                    rho, pv = stats.spearmanr(ok.score, ok[col])
                    cor.append({"programme": p, "cell_class": c, "donors": grp_name, "versus": col, "rho": rho,
                                "p": pv, "n": len(ok)})
    t = pd.DataFrame(rows)
    t["fdr"] = stats.false_discovery_control(t.p)
    t.to_csv(sp.P / "green2024_specificity.tsv", sep="\t", index=False)
    c = pd.DataFrame(cor)
    c.to_csv(sp.P / "green2024_stress_correlation.tsv", sep="\t", index=False)
    pd.set_option("display.width", 250)
    print(t.pivot_table(index="programme", columns="group", values="effect_sd").round(2)
          .reindex(columns=sp.GROUPS).to_string())
    print(c.pivot_table(index="programme", columns=["versus", "donors"], values="rho").round(2).to_string())


if __name__ == "__main__":
    {"pd": run_pd, "pd_nocore": lambda: run_pd(True), "psychad": run_psychad}[sys.argv[1]]()

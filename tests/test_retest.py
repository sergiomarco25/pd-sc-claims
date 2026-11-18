"""The engine's claim re-tests (claimstore.retest) and claim meta-analysis (claimstore.claim_meta) on a small synthetic
DE cache: a gene-set claim whose genes go up in cases is reproduced, the opposite claim contradicted, and a "no
difference" claim on unrelated genes reproduced; the stage contrast uses the stage of config/contrasts.toml.

Needs pandas and scipy (no anndata: gene-level tests read only the DE cache):  pytest tests/
"""
import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

SCRIPT = r'''
import json, sys
from collections import namedtuple
import numpy as np, pandas as pd
sys.path.insert(0, sys.argv[1])
from claimstore import retest, claim_meta
from claimstore import config as cfg

rng = np.random.default_rng(1)
ds, reg, cls = "synth", "R1", "Class A"
de = retest.DE / ds
(de / "de").mkdir(parents=True); (de / "input").mkdir()
n = 3000  # enough genes that 8 changed ones barely move library size
genes = [f"G{i}" for i in range(n)]
units = [f"d{i}|{reg}" for i in range(20)]
case = np.arange(20) < 10
stage = np.where(case, rng.integers(2, 7, 20), 0)
counts = pd.DataFrame(rng.poisson(50, (n, 20)), index=genes, columns=units)
counts.iloc[:8, case] *= 3                     # G0-G7 up in cases
counts.iloc[:8] = (counts.iloc[:8] * (1 + stage / 6)).round().astype(int)
counts.to_csv(de / "input" / f"{reg}__{retest.fname(cls)}.csv.gz")
pd.DataFrame({"unit": units, "region": reg, "disease": np.where(case, "case", "control"),
              "case_vs_control": np.where(case, "case", "ref"), retest.STAGE["columns"][0]: stage}).to_csv(
    de / "input" / "units.tsv", sep="\t", index=False)
lfc = np.r_[np.full(8, 1.5), rng.normal(0, 0.1, n - 8)]
p = np.r_[np.full(8, 1e-8), rng.uniform(0.05, 1, n - 8)]
for con in ("case_vs_control", retest.STAGE_CONTRAST):
    pd.DataFrame({"gene": genes, "logFC": lfc, "p": p, "padj": np.minimum(p * n, 1), "covariates": ""}).to_csv(
        de / "de" / f"{reg}__{retest.fname(cls)}__{con}.tsv.gz", sep="\t", index=False)
pd.DataFrame([{"file": "", "region": reg, "cell_class": cls, "contrast": "case_vs_control", "n_case": 10,
               "n_ref": 10, "n_genes": n, "n_fdr05": 8, "covariates": ""}]).to_csv(de / "de_log.tsv", sep="\t",
                                                                                    index=False)
Ctx = namedtuple("Ctx", "claim_id regions cell_classes direction contrast test_template cell_state null_claim")
sig = lambda gs: pd.DataFrame({"gene": gs, "kind": "named_genes", "signature_id": "s1", "variant": "both"})
out = {}
for name, gs, direction, con, null in [("up", genes[:8], "up", "case_vs_control", False),
                                       ("down", genes[:8], "down", "case_vs_control", False),
                                       ("null", genes[100:110], "up", "case_vs_control", True),
                                       ("stage", genes[:8], "associated_pos", retest.STAGE_CONTRAST, False),
                                       ("one", genes[:1], "up", "case_vs_control", False)]:
    c = Ctx("x-C01", reg, cls, direction, con, "de_pseudobulk", "", str(null))
    r = retest.test_claim(c, ds, sig(gs))[0]
    out[name] = {k: r.get(k) for k in ("outcome", "set_test", "template_used")}
out["no_genetics"] = retest.test_claim(Ctx("x-C02", reg, cls, "up", "case_vs_control", "genetic_enrichment", "", "False"),
                                       ds, sig([]), pd.Series({"statement": "s", "modality": "snRNA"}))[0]["outcome"]
g = pd.DataFrame({"outcome": ["reproduced"] * 3 + ["direction_only"], "z_oriented": [3.0, 2.5, 2.0, 1.0],
                  "n_units": [20, 30, 40, 10], "dataset": list("abcd")})
m = claim_meta.combine(g)
out["meta"] = [int(m["k"]), bool(m["p"] < 0.05), claim_meta.robustness(m)]
out["stage_label"] = retest.STAGE_LABEL
print(json.dumps(out))
'''


def run_synthetic():
    with tempfile.TemporaryDirectory() as tmp:
        env = {**os.environ}
        for k in list(env):
            if k.endswith(("_WORK", "_DATA_ROOT")):
                env.pop(k)
        sys.path.insert(0, str(ROOT / "src"))
        from claimstore import config as cfg
        env[f"{cfg.ENV_PREFIX}_WORK"] = tmp
        r = subprocess.run([sys.executable, "-I", "-c", SCRIPT, str(ROOT / "src")], capture_output=True, text=True,
                           env=env, cwd=tmp)
        assert r.returncode == 0, r.stderr[-3000:]
        return json.loads(r.stdout.strip().splitlines()[-1])


def test_gene_level_claims_on_a_synthetic_dataset():
    out = run_synthetic()
    assert out["up"]["outcome"] == "reproduced", out
    assert out["up"]["set_test"] == "donor signature score (Mann-Whitney)" and out["up"]["template_used"] == "de_genes"
    assert out["down"]["outcome"] == "contradicted", out
    assert out["null"]["outcome"] == "reproduced", out          # no difference claimed, none found
    assert out["stage"]["outcome"] == "reproduced", out
    assert out["stage"]["set_test"] == f"donor score ~ {out['stage_label']} (Spearman)"
    assert out["one"]["outcome"] == "reproduced" and out["one"]["set_test"] == "single gene (FDR)"
    assert out["no_genetics"] in ("not_testable", "reproduced", "contradicted", "direction_only", "not_reproduced")
    assert out["meta"] == [4, True, "robust"], out


if __name__ == "__main__":
    test_gene_level_claims_on_a_synthetic_dataset()
    print("PASS")

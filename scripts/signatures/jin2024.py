"""Signatures for jin2024 (Jin et al. 2024, Sci Adv; temporal cortex snRNA, 3 control / 5 LBD).

Reads (extracted from sciadv.adk3700_tables_s2_s3_s5_s6_s9_s10.zip into <supplement>/extracted/):
  adk3700_table_s5.xlsx  table S5: post-mortem brain DEGs per cell type, LBD vs control (MAST)
  adk3700_table_s6.xlsx  table S6: IPA canonical pathways + toxicity lists per cell type (Molecules)
(tables S2/S3 are organoid DEGs/pathways and are not used.)
Writes claims/signatures/jin2024.tsv and claims/signatures/jin2024_unresolved.tsv.
Run: ~/miniconda3/envs/rapids_singlecell/bin/python -I scripts/signatures/jin2024.py

Notes
- Table S5 is already the significant set (all FDR < 0.05). EX has 188 up / 966 down, matching the
  text (claim C03); all named genes of C04-C08 are present with the stated direction.
- IPA 'Molecules' are DEG members of each pathway only (IPA gene sets are proprietary). The suffix
  ' (includes others)' that IPA appends to some symbols is removed.
"""
from pathlib import Path
import zipfile
import pandas as pd
import yaml

PAPER = "jin2024"
ROOT = Path(__file__).resolve().parents[2]

import sys  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

SUP = cfg.RAW / "supplements" / PAPER
EXT = SUP / "extracted"
ZIP = SUP / "sciadv.adk3700_tables_s2_s3_s5_s6_s9_s10.zip"
OUT = ROOT / "claims" / "signatures"
COLS = ["signature_id", "paper", "claim_ids", "kind", "cell_type", "cell_state", "comparison",
        "gene", "rank", "stat", "padj", "source", "selection"]
C = lambda *n: [f"{PAPER}-C{int(i):02d}" for i in n]  # noqa: E731
CT = {"EX": "Excitatory neurons", "IN": "Inhibitory neurons", "OLG": "Oligodendrocytes", "OPC": "OPCs",
      "AS": "Astrocytes", "MG": "Microglia", "EC": "Endothelial"}
COMP = "LBD vs control"

if not (EXT / "adk3700_table_s5.xlsx").exists():  # extract the zip into the paper's own folder
    with zipfile.ZipFile(ZIP) as z:
        for m in z.namelist():
            if m.endswith(".xlsx") and "/" not in m and ".." not in m:
                z.extract(m, EXT)
X5, X6 = EXT / "adk3700_table_s5.xlsx", EXT / "adk3700_table_s6.xlsx"

rows, unresolved, DATE_GENES = [], [], set()


def gname(g):
    if hasattr(g, "strftime"):  # Excel date-corrupted symbol: keep as Excel displays it, flag
        s = f"{g.day}-{g.strftime('%b')}"
        DATE_GENES.add(s)
        return s
    return str(g).replace(" (includes others)", "").strip()


def add(name, claims, kind, ct, state, comp, genes, stats, padjs, source, selection):
    genes = [gname(g) for g in genes]
    for i, g in enumerate(genes):
        rows.append(dict(signature_id=f"{PAPER}__{name}", paper=PAPER, claim_ids=";".join(claims),
                         kind=kind, cell_type=ct, cell_state=state, comparison=comp, gene=g, rank=i + 1,
                         stat="" if stats is None else stats[i], padj="" if padjs is None else padjs[i],
                         source=source, selection=selection))


def unres(claim, term, reason):
    unresolved.append(dict(claim_id=claim, term=term, reason=reason))


# ---------------- table S5: DEGs ----------------
s5 = {}
for sheet, d in pd.read_excel(X5, sheet_name=None).items():
    d.columns = [str(c).strip() for c in d.columns]
    ctk = sheet.split()[-1]
    d = d[d["Gene Name"].notna() & (d.FDR < 0.05)].copy()
    s5[ctk] = d
    cl = C(3, 4, 5, 6, 7, 8, 9, 10, 11) if ctk == "EX" else C(11)
    for kind, sub in (("de_up", d[d.LogFC > 0]), ("de_down", d[d.LogFC < 0])):
        sub = sub.reindex(sub.LogFC.abs().sort_values(ascending=False).index)
        if len(sub):
            add(f"{ctk}_LBD_vs_control_{kind.split('_')[1]}", cl, kind, CT[ctk], "", COMP, sub["Gene Name"],
                sub.LogFC.round(4).tolist(), sub.FDR.tolist(), f"{X5.name}|{sheet}",
                "table S5 (MAST, sample random effect): FDR<0.05 (authors' |log2FC|>=0.1 already applied); "
                "ranked by |LogFC|")


# ---------------- table S6: IPA pathways ----------------
def ipa(sheet):
    raw = pd.read_excel(X6, sheet_name=sheet, header=None)
    out, section, header = [], None, None
    for _, r in raw.iterrows():
        first = str(r.iloc[0]).strip()
        if first in ("Ingenuity Canonical Pathways", "Ingenuity Toxicity Lists"):
            section, header = first, [str(x).strip() for x in r.values]
            continue
        if section is None or pd.isna(r.iloc[0]):
            continue
        rec = dict(zip(header, r.values))
        out.append(dict(section=section, pathway=first, mlog10p=rec["-log(p-value)"],
                        genes=[g.strip() for g in str(rec["Molecules"]).split(",") if g.strip() and g != "nan"]))
    return pd.DataFrame(out)


S6 = {sh.split()[-1]: ipa(sh) for sh in pd.ExcelFile(X6).sheet_names}

# C09: canonical pathways named for EX
lfc_ex = s5["EX"].set_index("Gene Name").LogFC
c09 = ["Calcium Signaling", "Synaptic Long Term Potentiation", "Synaptogenesis Signaling Pathway",
       "SNARE Signaling Pathway", "Clathrin-mediated Endocytosis Signaling", "RHOGDI Signaling",
       "Protein Kinase A Signaling", "α-Adrenergic Signaling", "Unfolded protein response",
       "Glycogen Degradation II", "Protein Ubiquitination Pathway"]
ex = S6["EX"]
for pw in c09:
    r = ex[(ex.section == "Ingenuity Canonical Pathways") & (ex.pathway == pw)]
    assert len(r) == 1, pw
    g = pd.Series([gname(x) for x in r.genes.iloc[0]]).drop_duplicates()
    g = g.reindex(g.map(lambda x: -abs(lfc_ex.get(x, 0))).sort_values(kind="stable").index)
    short = "".join(ch if ch.isalnum() else "_" for ch in pw.replace("α", "alpha")).strip("_")
    add(f"EX_IPA_{short}", C(9), "pathway_genes", CT["EX"], "", COMP, g,
        [round(lfc_ex[x], 4) if x in lfc_ex else "" for x in g], None, f"{X6.name}|Table S6 EX",
        f"IPA canonical pathway '{pw}' (-log10 p {r.mlog10p.iloc[0]}): DEG members listed by IPA; "
        "ranked by |LogFC| in table S5 EX (stat = LogFC)")
    unres(C(9)[0], f"IPA: {pw}", "IPA gene set proprietary; signature holds only the EX DEG members; "
          "use closest GO/Reactome set for a full gene-set test")

# C10: toxicity lists for EX (mitochondrial dysfunction, lipid metabolism, oxidative stress)
tox = ex[(ex.section == "Ingenuity Toxicity Lists")
         & ex.pathway.str.contains("Mitochondri|Fatty Acid|Lipid|Oxidative", case=False)]
for _, r in tox.iterrows():
    g = pd.Series([gname(x) for x in r.genes]).drop_duplicates()
    g = g.reindex(g.map(lambda x: -abs(lfc_ex.get(x, 0))).sort_values(kind="stable").index)
    short = "".join(ch if ch.isalnum() else "_" for ch in r.pathway)[:40].strip("_")
    add(f"EX_IPAtox_{short}", C(10), "pathway_genes", CT["EX"], "", COMP, g,
        [round(lfc_ex[x], 4) if x in lfc_ex else "" for x in g], None, f"{X6.name}|Table S6 EX",
        f"IPA toxicity list '{r.pathway}' (-log10 p {r.mlog10p}): DEG members; ranked by |LogFC| in table S5 EX")
unres(C(10)[0], "IPA toxicity lists (Fig. 3G): mitochondrial dysfunction, lipid metabolism, oxidative stress",
      "IPA lists proprietary; signatures hold EX DEG members of matching toxicity lists in table S6 only")

# C11: pathways perturbed across cell types -> union over cell types
c11 = ["RHOGDI Signaling", "Protein Kinase A Signaling", "Synaptogenesis Signaling Pathway",
       "Signaling by Rho Family GTPases", "Estrogen Receptor Signaling", "Protein Ubiquitination Pathway"]
for pw in c11:
    hits = {ct: t[(t.section == "Ingenuity Canonical Pathways") & (t.pathway == pw)] for ct, t in S6.items()}
    hits = {ct: h for ct, h in hits.items() if len(h)}
    cnt = pd.Series([gname(x) for h in hits.values() for x in set(h.genes.iloc[0])]).value_counts()
    cnt = cnt.reindex(sorted(cnt.index, key=lambda x: (-cnt[x], x)))
    short = "".join(ch if ch.isalnum() else "_" for ch in pw).strip("_")
    add(f"All_IPA_{short}", C(11), "pathway_genes", "All", "", COMP, cnt.index, cnt.tolist(), None,
        f"{X6.name}|Table S6 " + ",".join(hits),
        f"IPA canonical pathway '{pw}': union of DEG members over cell types {','.join(hits)}; "
        "stat = number of cell types listing the gene")
    unres(C(11)[0], f"IPA: {pw}", "IPA gene set proprietary; signature holds only DEG members; use GO/Reactome equivalent")

# ---------------- named genes ----------------
named = [
    ("EX_SNCA", C(2), "Excitatory neurons", "", "cell-type expression", ["SNCA"]),
    ("EX_synaptic_up_named", C(4), "Excitatory neurons", "", COMP, ["GRIA4", "PSD3", "CACNB2", "KCNMB2", "GRIN2B"]),
    ("EX_synaptic_down_named", C(5), "Excitatory neurons", "", COMP, ["PCLO", "NCALD", "ATP2B2"]),
    ("EX_mito_down_named", C(6), "Excitatory neurons", "", COMP, ["PRKN", "ABCB10", "MFN1", "COX15"]),
    ("EX_lipid_down_named", C(7), "Excitatory neurons", "", COMP, ["FADS1", "PEX5L"]),
    ("EX_oxstress_named", C(8), "Excitatory neurons", "", COMP + " (OXR1, IDH2 down; PRKAR2B up)",
     ["OXR1", "IDH2", "PRKAR2B"]),
]
for name, cl, ct, st, comp, genes in named:
    add(name, cl, "named_genes", ct, st, comp, genes, None, None, "text", "text list")

# ---------------- write + checks ----------------
OUT.mkdir(parents=True, exist_ok=True)
sig = pd.DataFrame(rows, columns=COLS)
un = pd.DataFrame(unresolved, columns=["claim_id", "term", "reason"])
assert sig.gene.str.len().gt(0).all() and not sig.gene.isin(["nan", "None"]).any()
assert (sig.groupby("signature_id")["rank"].min() == 1).all()
dup = sig.duplicated(["signature_id", "gene"])
if dup.any():
    print("WARNING duplicated genes within signature (kept):", sig[dup][["signature_id", "gene"]].values.tolist())
sig.to_csv(OUT / f"{PAPER}.tsv", sep="\t", index=False)
un.to_csv(OUT / f"{PAPER}_unresolved.tsv", sep="\t", index=False)

print(sig.groupby(["kind", "signature_id"]).size().to_string())
print("\nsignatures by kind:", sig.drop_duplicates("signature_id").kind.value_counts().to_dict())
susp = sig[sig.gene.str.match(r"^\d+-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$|^(Jan|Feb|Mar|Apr|Sep|Oct|Nov|Dec)-\d+$|^[\d.]+$")]
print("possible Excel date-corrupted / numeric symbols:", sorted(set(susp.gene) | DATE_GENES))
print("multi-symbol IPA entries (kept as given):", sorted(sig[sig.gene.str.contains("/")].gene.unique()))

claims = yaml.safe_load(open(ROOT / "claims" / f"{PAPER}.yaml"))
need = {c["id"] for c in claims if c["test_template"] in ("state_signature", "geneset_shift", "de_pseudobulk")
        or (c["test_template"].startswith("composition") and (c.get("cell_state") or "").strip())}
in_sig = {x for s in sig.claim_ids.unique() for x in s.split(";")}
in_un = set(un.claim_id)
print(f"\ncoverage: {len(need)} claims need signatures; {len(need & in_sig)} in signatures; "
      f"{len(need & in_un - in_sig)} unresolved only; missing: {sorted(need - in_sig - in_un)}")

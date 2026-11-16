"""Consensus statements: published claims that make the same point, grouped and judged together.

Usage: python -I scripts/42_consensus.py
Input: claims/consensus_statements.tsv (curated: each statement lists the claims that assert it and those that assert
the opposite, a summary of the evidence from the joint analyses and a verdict) and claims/claim_evidence.tsv (tier of
every claim). Writes claims/consensus_summary.tsv, manuscript/table1_consensus.md and figures/B16_consensus.
"""
from pathlib import Path
import sys
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))
from pdrep.plot.style import apply_style, save_figure, ROBUSTNESS_COLORS, INK, MUTED, GRID  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

TIERS = ["1 established", "2 supported", "3 unresolved", "4 contradicted", "5 untested"]
TCOL = {"1 established": ROBUSTNESS_COLORS["robust"], "2 supported": ROBUSTNESS_COLORS["supported"],
        "3 unresolved": ROBUSTNESS_COLORS["not_replicated"], "4 contradicted": ROBUSTNESS_COLORS["contradicted"],
        "5 untested": "#f2f2f2"}
TLAB = {"1 established": "established", "2 supported": "supported", "3 unresolved": "unresolved",
        "4 contradicted": "contradicted", "5 untested": "untested"}
VCOL = {"established": "#313695", "supported": "#5b8fc7", "supported when pooled": "#5b8fc7",
        "excitatory loss supported": "#5b8fc7", "partly supported": "#9fb7d6",
        "AGTR1 loss yes, selectivity no": "#9fb7d6", "not supported as proportions": "#f4a259",
        "not supported": "#f4a259", "largely not supported": "#f4a259", "contradicted": "#c0282d"}
PAPER = {"smajic2022": "Smajić", "kamath2022": "Kamath", "lee2023": "Lee 2023", "martirosyan2024": "Martirosyan",
         "wang2024": "Wang", "adams2024": "Adams", "ma2025": "Ma", "volpato2025": "Volpato", "menon2026": "Menon",
         "barbareyes2026": "Barba-Reyes 2026", "barbareyes2025": "Barba-Reyes 2025", "xu2023": "Xu",
         "prashant2024": "Prashant", "feleke2021": "Feleke", "zhu2024": "Zhu", "chia2025": "Chia", "jin2024": "Jin",
         "glock2026": "Glock", "lee2026": "Lee 2026", "clarence2025": "Clarence", "dehestani2024": "Dehestani"}


def split(s):
    return [x for x in str(s).split(";") if x and x != "nan"]


def main():
    cs = pd.read_csv(ROOT / "claims" / "consensus_statements.tsv", sep="\t")
    ev = pd.read_csv(ROOT / "claims" / "claim_evidence.tsv", sep="\t").set_index("claim_id")
    rows = []
    for r in cs.itertuples():
        sup, opp = split(r.supporting_claims), split(r.opposing_claims)
        missing = [c for c in sup + opp if c not in ev.index]
        assert not missing, (r.id, missing)
        row = {"id": r.id, "statement": r.statement, "short": r.short, "verdict": r.verdict, "our_evidence": r.our_evidence,
               "n_claims": len(sup) + len(opp), "n_supporting": len(sup), "n_opposing": len(opp),
               "papers": ", ".join(sorted({PAPER.get(ev.loc[c, "paper"], ev.loc[c, "paper"]) for c in sup + opp}))}
        row["n_papers"] = len(row["papers"].split(", "))
        for side, ids in [("sup", sup), ("opp", opp)]:
            t = ev.loc[ids, "tier"].value_counts() if ids else pd.Series(dtype=int)
            for k in TIERS:
                row[f"{side}_{TLAB[k]}"] = int(t.get(k, 0))
        rows.append(row)
    s = pd.DataFrame(rows)
    s.to_csv(ROOT / "claims" / "consensus_summary.tsv", sep="\t", index=False)
    # manuscript table
    lines = ["| | Consensus statement | Claims (papers) | Supporting claims held / tested | Evidence from the joint analysis | Verdict |",
             "|---|---|---|---|---|---|"]
    for r in s.itertuples():
        held = r.sup_established + r.sup_supported
        tested = r.n_supporting - r.sup_untested
        lines.append(f"| {r.id[2:]} | {r.statement} | {r.n_claims} ({r.n_papers}) | {held}/{tested} | {r.our_evidence} | **{r.verdict}** |")
    (ROOT / "manuscript" / "table1_consensus.md").write_text("\n".join(lines) + "\n")
    fig_(s)
    print(s[["id", "short", "n_claims", "n_papers", "sup_established", "sup_supported", "sup_unresolved",
             "sup_contradicted", "sup_untested", "verdict"]].to_string(index=False))


def fig_(s):
    apply_style(font__size=7, axes__linewidth=0.6, axes__edgecolor=INK)
    fig, ax = plt.subplots(figsize=(7.2, 0.3 * len(s) + 1.2))
    y = np.arange(len(s))
    mx = (s.n_supporting).max()
    for i, r in enumerate(s.itertuples()):
        x = 0
        for k in TIERS:
            n = getattr(r, f"sup_{TLAB[k]}")
            if n:
                ax.add_patch(Rectangle((x, i - 0.35), n, 0.7, color=TCOL[k], lw=0))
                x += n
        x = 0  # opposing claims drawn to the left
        for k in TIERS:
            n = getattr(r, f"opp_{TLAB[k]}")
            if n:
                ax.add_patch(Rectangle((-x - n, i - 0.35), n, 0.7, color=TCOL[k], lw=0, hatch="////", ec="white"))
                x += n
        ax.text(mx + 1, i, r.verdict, va="center", fontsize=6.2, color=VCOL.get(r.verdict, INK), fontweight="bold")
        ax.text(-7.5, i, f"{r.short}  ({r.n_papers} paper{'s' if r.n_papers > 1 else ''})", va="center", ha="right", fontsize=6.4)
    ax.axvline(0, color=INK, lw=0.6)
    ax.set_xlim(-7, mx + 14)
    ax.set_ylim(len(s) - 0.5, -0.8)
    ax.set_yticks([])
    ax.set_xticks([-5, 0, 5, 10, 15, 20, 25, 30], ["5", "0", "5", "10", "15", "20", "25", "30"])
    ax.set_xlabel("claims asserting the opposite (left, hatched) and the statement (right), by cross-dataset evidence")
    for sp_ in ["top", "right", "left"]:
        ax.spines[sp_].set_visible(False)
    ax.legend(handles=[Rectangle((0, 0), 1, 1, color=TCOL[k], label=TLAB[k]) for k in TIERS], frameon=False,
              fontsize=6, ncol=5, loc="lower center", bbox_to_anchor=(0.45, 1.0))
    ax.text(mx + 1, -0.75, "verdict", fontsize=6.2, color=MUTED)
    fig.subplots_adjust(left=0.36, right=0.98, top=0.9, bottom=0.12)
    save_figure(fig, ROOT / "figures" / "B16_consensus")
    plt.close(fig)


if __name__ == "__main__":
    main()

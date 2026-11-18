"""Claim gene signatures (claimstore.signatures): named genes from the claims, replaced by curated signatures."""
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore import signatures as S  # noqa: E402

CLAIMS = [
    {"id": "p1-C01", "paper": "p1", "genes": ["GENE1", "GENE2", "GENE1"], "cell_type": ["Microglia"]},
    {"id": "p1-C02", "paper": "p1", "genes": [], "cell_type": "All"},
    {"id": "p2-C01", "paper": "p2", "genes": "GENE3; GENE4", "cell_type": "All"},
]


def test_named_genes_from_claims():
    sig = S.build(CLAIMS, pd.DataFrame(columns=S.COLUMNS))
    assert list(sig.columns) == S.COLUMNS
    assert sig.groupby("claim_ids").gene.apply(list).to_dict() == {"p1-C01": ["GENE1", "GENE2"],
                                                                   "p2-C01": ["GENE3", "GENE4"]}
    assert set(sig.kind) == {"named_genes"} and sig.cell_type[sig.claim_ids.eq("p1-C01")].eq("Microglia").all()
    assert sig.selection.eq("named in the claim text").all()


def test_curated_signature_replaces_named_genes(tmp_path):
    (tmp_path / "p2.tsv").write_text("signature_id\tpaper\tclaim_ids\tkind\tcell_type\tgene\n"
                                     "p2_up\tp2\tp2-C01;p2-C09\tde_up\tAll\tGENE9\n")
    (tmp_path / "p2_unresolved.tsv").write_text("claim_id\treason\np2-C05\tno table\n")
    cur = S.curated(tmp_path)
    assert list(cur.columns) == S.COLUMNS and cur.gene_reported.tolist() == ["GENE9"]
    sig = S.build(CLAIMS, cur)
    assert sig[sig.claim_ids.str.contains("p2-C01")].gene.tolist() == ["GENE9"]
    assert sig[sig.claim_ids.eq("p1-C01")].gene.tolist() == ["GENE1", "GENE2"]
    assert S.curated(tmp_path / "missing").empty

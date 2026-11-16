"""volpato2025 — Volpato et al. 2025 medRxiv (GSE329625), full-length single-nucleus RNA-seq (ICELL8 / SMART-seq,
Takara Mappa pipeline, STAR + Subread) of ventral SNpc.

Sources (raw/volpato2025/):
  GSE329625_processed_gene_matrix.txt.gz: genes x nuclei, tab-separated, header = nucleus IDs prefixed "X"
      (R make.names). Integer read counts (Subread): exon-derived counts for all protein-coding genes plus
      intron-derived counts for the authors' top 10,000 variable protein-coding genes, summed per gene (GEO
      data_processing). Read counts, not UMIs and not TPM.
  GSE329625_full_metadata.csv.gz: per-nucleus Sample_v2 (donor/library), Disease, Sex, CellType, CellSubType.
Disease and sex are taken from full_metadata. The GEO soft file pairs library names with disease/sex/age
differently (e.g. library PD703 listed as Control), so its per-sample ages are not used (age = NaN).
"""
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from pdrep.harmonise import RAW

D = RAW / "volpato2025"
AUTHOR_TO_CLASS = {
    **{f"DaN_{i}": "DA neurons" for i in range(4)},
    **{f"ODC_{i}": "Oligodendrocytes" for i in range(3)},
    **{f"Astrocyte_{i}": "Astrocytes" for i in range(4)},
    **{f"Microglia_{i}": "Microglia" for i in range(3)},
    "OPC": "OPCs", "GABA": "Inhibitory neurons", "Tcells": "T cells",
}
DISEASE = {"CTR": "control", "PD_B5-6": "PD", "PD_B3-4": "PD", "ILBD_B3-4": "iLBD", "ILBD_B1-2": "iLBD"}


def load():
    m = pd.read_csv(D / "GSE329625_full_metadata.csv.gz", index_col=0, dtype={"Sample_v2": str})
    mats, genes = [], []
    with pd.read_csv(D / "GSE329625_processed_gene_matrix.txt.gz", sep="\t", index_col=0, chunksize=2000) as it:
        for ch in it:
            v = ch.to_numpy()
            assert np.array_equal(v, np.round(v)), "non-integer counts"
            mats.append(sp.csr_matrix(v.astype(np.float32)))
            genes.extend(ch.index.astype(str))
            cells = [c[1:] if c.startswith("X") else c for c in ch.columns.astype(str)]
    X = sp.vstack(mats).T.tocsr()
    missing = set(cells) - set(m.index)
    assert not missing, f"{len(missing)} matrix nuclei not in metadata, e.g. {sorted(missing)[:3]}"
    m = m.loc[cells]
    obs = pd.DataFrame({
        "donor": m.Sample_v2.values, "sample": m.Sample_v2.values, "region": "SN",
        "disease": m.Disease.map(DISEASE).values, "disease_raw": m.Disease.values,
        "braak_lb_group": m.Disease.str.extract(r"B(\d-\d)")[0].values,
        "sex": m.Sex.values, "age": np.nan, "pmi": np.nan,
        # stage group midpoints (B1-2 -> 1.5, B3-4 -> 3.5, B5-6 -> 5.5); controls left NaN (stage 0 in contrasts)
        "braak_lb": m.Disease.str.extract(r"B(\d)-(\d)").astype(float).mean(axis=1).values,
        "author_cell_type": m.CellSubType.values, "author_cell_type_broad": m.CellType.values,
        "percent_mt_author": m["percent.mt"].values, "platform": "ICELL8 full-length (read counts)",
    }, index=m.index)
    assert obs.disease.notna().all()
    return ad.AnnData(X, obs=obs, var=pd.DataFrame(index=genes))

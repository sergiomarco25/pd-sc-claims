"""Is the heat-shock response an active transcriptional programme? HSF1-motif accessibility in paired multiome data.

Usage: python -I scripts/44_hsf1_accessibility.py
Data: Adams et al. 2024 midbrain multiome (Cell Ranger ARC matrices, one library per donor; RNA and ATAC from the same
nuclei). Nuclei are labelled with the harmonised cell classes (harmonised/adams2024.h5ad, index <sample>_<barcode>).
For every library, ATAC peak counts are summed per cell class. Peaks are scanned in hg38 for the three-unit heat-shock
element (nGAAnnTTCnnGAAn, either orientation: GAA..TTC..GAA or TTC..GAA..TTC) and for a composition-matched control
k-mer (AAG..CTT..AAG or CTT..AAG..CTT). Per donor x class, HSF1 accessibility = log of the share of ATAC counts in
HSE peaks over the share in control-k-mer peaks. Compared between PD and aged controls (young controls excluded) and
correlated across donors with the RNA heat-response score of the same donor and class
(results/biology/donor_scores.tsv.gz). Writes results/biology/hsf1_accessibility.tsv and hsf1_accessibility_tests.tsv.
"""
from pathlib import Path
import re
import sys
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

RAW = cfg.RAW / "adams2024"
HARM = cfg.HARMONISED / "adams2024.h5ad"
B = cfg.RESULTS / "biology"
FA = cfg.GENOME_FA
HSE = re.compile(r"(?=(GAA..TTC..GAA|TTC..GAA..TTC))")
CTRL = re.compile(r"(?=(AAG..CTT..AAG|CTT..AAG..CTT))")
CLASSES = ["Oligodendrocytes", "Astrocytes", "Microglia", "OPCs"]


class Fasta:
    def __init__(self, path):
        self.f = open(path, "rb")
        self.idx = {r[0]: (int(r[1]), int(r[2]), int(r[3]), int(r[4]))
                    for r in (l.split("\t") for l in open(str(path) + ".fai"))}

    def fetch(self, chrom, start, end):
        n, off, lb, lw = self.idx[chrom]
        end = min(end, n)
        a = off + (start // lb) * lw + start % lb
        b = off + (end // lb) * lw + end % lb
        self.f.seek(a)
        return self.f.read(b - a).replace(b"\n", b"").decode().upper()


def peak_flags(peaks, fa, cache):
    hse, ctl = np.zeros(len(peaks), bool), np.zeros(len(peaks), bool)
    for i, p in enumerate(peaks):
        if p in cache:
            hse[i], ctl[i] = cache[p]
            continue
        chrom, rng = p.split(":")
        s, e = (int(x) for x in rng.split("-"))
        if chrom not in fa.idx:
            cache[p] = (False, False)
            continue
        seq = fa.fetch(chrom, s, e)
        cache[p] = (bool(HSE.search(seq)), bool(CTRL.search(seq)))
        hse[i], ctl[i] = cache[p]
    return hse, ctl


def main():
    obs = ad.read_h5ad(HARM, backed="r").obs[["donor", "disease", "disease_raw", "age_group", "cell_class"]].copy()
    fa = Fasta(FA)
    cache, rows = {}, []
    for f in sorted(RAW.glob("*_filtered_feature_bc_matrix.h5")):
        sample = f.name.split("_")[1]
        with h5py.File(f, "r") as h:
            m = h["matrix"]
            ft = m["features/feature_type"][:]
            ids = m["features/id"][:].astype(str)
            bc = m["barcodes"][:].astype(str)
            X = sp.csc_matrix((m["data"][:], m["indices"][:], m["indptr"][:]), shape=m["shape"][:]).tocsr()
        pk = ft == b"Peaks"
        Xp = X[pk]
        peaks = ids[pk]
        hse, ctl = peak_flags(peaks, fa, cache)
        lab = obs.reindex([f"{sample}_{b}" for b in bc])
        for c in CLASSES:
            cols = np.where(lab.cell_class.astype(str).values == c)[0]
            if len(cols) < 20:
                continue
            v = np.asarray(Xp[:, cols].sum(1)).ravel()
            tot = v.sum()
            rows.append({"sample": sample, "cell_class": c, "nuclei": len(cols),
                         "donor": lab.donor.iloc[cols[0]], "disease": lab.disease.iloc[cols[0]],
                         "age_group": lab.age_group.iloc[cols[0]], "atac_counts": tot,
                         "hse_share": v[hse].sum() / tot, "ctrl_share": v[ctl].sum() / tot,
                         "hse_peaks": int(hse.sum()), "ctrl_peaks": int(ctl.sum()), "peaks": len(peaks)})
        print(sample, len(peaks), int(hse.sum()), int(ctl.sum()), flush=True)
    d = pd.DataFrame(rows)
    d["hsf1_access"] = np.log(d.hse_share / d.ctrl_share)
    s = pd.read_csv(B / "donor_scores.tsv.gz", sep="\t", low_memory=False)
    s = s[(s.dataset == "adams2024") & (s.programme == "heat response")][["donor", "cell_class", "score"]]
    d = d.merge(s, on=["donor", "cell_class"], how="left")
    d.to_csv(B / "hsf1_accessibility.tsv", sep="\t", index=False)
    tests = []
    for c, x in d.groupby("cell_class"):
        aged = x[~x.age_group.astype(str).str.lower().str.contains("young")]
        a, r = aged[aged.disease.eq("PD")].hsf1_access, aged[aged.disease.eq("control")].hsf1_access
        u = stats.mannwhitneyu(a, r) if len(a) >= 3 and len(r) >= 3 else None
        y = x.dropna(subset=["score"])
        rho = stats.spearmanr(y.hsf1_access, y.score) if len(y) >= 8 else None
        tests.append({"cell_class": c, "n_pd": len(a), "n_aged_ctrl": len(r), "pd_minus_ctrl": a.mean() - r.mean(),
                      "p_pd_vs_ctrl": u.pvalue if u else np.nan, "n_corr": len(y),
                      "rho_access_vs_rna_heat": rho[0] if rho else np.nan, "p_rho": rho[1] if rho else np.nan})
    t = pd.DataFrame(tests)
    t.to_csv(B / "hsf1_accessibility_tests.tsv", sep="\t", index=False)
    print(t.round(4).to_string(index=False))


if __name__ == "__main__":
    main()

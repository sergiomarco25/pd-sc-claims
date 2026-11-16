"""feleke2021 — Feleke et al. 2021 Acta Neuropathol, GSE178146 (10x v3 snRNA-seq, anterior cingulate cortex,
Parkinson's UK Tissue Bank; 7 control, 7 PD, 7 PDD, 7 DLB; one library per donor).

Sources:
  raw/feleke2021/GSE178146_RAW.tar -> extracted/<GSM>_<id>/<id>_raw_feature_bc_matrix/  Cell Ranger 3.0.2 *raw*
      (unfiltered) matrices, GRCh38 / Ensembl 93. No filtered matrices or cell barcodes are public.
  supplements/feleke2021/401_2021_2343_MOESM2_ESM.xlsx  "Sample information" (disease group, sex, age at onset/death,
      disease duration, PMI, Braak aSN stage, Braak tau, Thal, ACG aSN score, Lewy pathology category) and
      "Single-nucleus RNA-seq" (authors' viable nuclei per sample, used only for comparison).
Cell calling (re-implemented here, as the paper did with DropletUtils 1.6.1 emptyDrops, lower = 300 UMI):
  * ambient profile = summed counts of barcodes with total UMI <= 300, Good-Turing-like smoothing (genes unseen in the
    ambient pool share the mass of singleton genes); Dirichlet-multinomial concentration alpha fitted by maximum
    likelihood on those ambient barcodes;
  * barcodes with total > 300 get the DM log-likelihood under the ambient profile; Monte Carlo p-values from 10,000 (DropletUtils default)
    simulated ambient droplets (Polya-urn sampling, compared at identical totals); BH FDR <= 0.001;
  * barcodes at or above the barcodeRanks *inflection* point are retained without testing (emptyDrops `retain`
    argument; DropletUtils' default is the knee). Reason: in these snRNA libraries the knee lands at 2,000-12,000 UMI
    (retaining only 100-2,000 barcodes), and in some libraries (e.g. PD294, PD332) the <=300-UMI ambient pool is so
    nucleus-like that barcodes with ~1,000 UMI are *more* probable under the ambient model than simulated ambient
    droplets, so emptyDrops calls none of them (sig = 0) and the sample would almost vanish. With retain = knee the
    first build (job 41861760) recovered 28-88 nuclei for such samples vs 3,000-4,500 in the paper.
  Then the uniform QC used in all loaders: >= 200 genes, < 10% mito UMIs (the paper used < 5% mito and DoubletFinder;
  doublets are not removed here). Nuclei per sample vs the paper are printed in the log.
Donor IDs: GEO uses C36/C48, the supplement C036/C048 -> donor uses the supplement form.
disease: Control / PD / PDD / DLB as given (PDD and DLB are the paper's clinical groups). braak_lb = Braak aSN stage.
No public author labels -> marker clustering in finalise().
"""
import re
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp
from scipy.special import gammaln
from scipy.optimize import minimize_scalar
from scipy.interpolate import make_smoothing_spline
import anndata as ad
from pdrep.harmonise import RAW

AUTHOR_TO_CLASS = {}
D = RAW / "feleke2021" / "extracted"
SUPP = RAW / "supplements" / "feleke2021" / "401_2021_2343_MOESM2_ESM.xlsx"
LOWER, FDR, N_SIM, SEED = 300, 0.001, 10000, 0  # niters = DropletUtils default (2,000 caps p at 1/2001: BH then fails)


def knee_point(tot, lower=LOWER, exclude_from=50):
    """DropletUtils::barcodeRanks: smoothing spline of log10(total) on log10(mid-rank) over unique totals > lower
    (top 50 points excluded); inflection = steepest point, knee = max curvature between the flattest point to its left
    and the inflection. Returns (knee, inflection) in UMI."""
    t = np.sort(tot[tot > 0])
    vals, cnt = np.unique(t, return_counts=True)
    vals, cnt = vals[::-1], cnt[::-1]
    mid = np.cumsum(cnt) - (cnt - 1) / 2
    keep = vals > lower
    x, y = np.log10(mid[keep])[exclude_from:], np.log10(vals[keep])[exclude_from:]
    spl = make_smoothing_spline(x, y)
    d1, d2 = spl.derivative(1)(x), spl.derivative(2)(x)
    right = int(np.argmin(d1))
    left = int(np.argmax(d1[:right + 1]))
    curv = d2 / (1 + d1 ** 2) ** 1.5
    return float(10 ** spl(x[left + np.argmin(curv[left:right + 1])])), float(10 ** spl(x[right]))


def ambient_profile(Xa):
    s = np.asarray(Xa.sum(0)).ravel().astype(float)
    n1 = (s == 1).sum()
    zero = s == 0
    p = s / s.sum()
    if zero.any():
        mass = max(n1, 1) / s.sum()
        p = p * (1 - mass)
        p[zero] = mass / zero.sum()
    return p


def dm_loglik(X, p, alpha):
    """Sum over nonzero genes of lgamma(y + a p) - lgamma(a p) - lgamma(y + 1) (terms depending only on T dropped)."""
    X = X.tocsr()
    ap = alpha * p[X.indices]
    v = gammaln(X.data + ap) - gammaln(ap) - gammaln(X.data + 1)
    return np.add.reduceat(v, X.indptr[:-1]) * (np.diff(X.indptr) > 0)


def fit_alpha(Xa, p):
    Xa = Xa[np.asarray(Xa.sum(1)).ravel() > 0]
    T = np.asarray(Xa.sum(1)).ravel()

    def nll(la):
        a = np.exp(la)
        return -(dm_loglik(Xa, p, a).sum() + (gammaln(a) - gammaln(T + a)).sum())
    return float(np.exp(minimize_scalar(nll, bounds=(-2, 15), method="bounded").x))


def simulate(p, alpha, tmax, totals_needed, n_sim=N_SIM, seed=SEED):
    """Polya-urn draws of ambient droplets; returns {T: sorted simulated logliks} for T in totals_needed."""
    rng = np.random.default_rng(seed)
    cp = np.cumsum(p)
    hist = np.zeros((n_sim, tmax), dtype=np.int32)
    counts = np.zeros((n_sim, len(p)), dtype=np.int32)
    ll = np.zeros(n_sim)
    out, need = {}, set(int(t) for t in totals_needed)
    rows = np.arange(n_sim)
    for t in range(tmax):
        fresh = rng.random(n_sim) < alpha / (t + alpha)
        g = np.empty(n_sim, dtype=np.int64)
        nf = fresh.sum()
        g[fresh] = np.minimum(np.searchsorted(cp, rng.random(nf) * cp[-1]), len(p) - 1)
        if t > 0 and nf < n_sim:
            j = rng.integers(0, t, n_sim - nf)
            g[~fresh] = hist[rows[~fresh], j]
        c = counts[rows, g]
        ll += np.log(c + alpha * p[g]) - np.log(c + 1)
        counts[rows, g] = c + 1
        hist[:, t] = g
        if t + 1 in need:
            out[t + 1] = np.sort(ll.copy())
    return out


def empty_drops(X, lower=LOWER):
    tot = np.asarray(X.sum(1)).ravel()
    amb = tot <= lower
    p = ambient_profile(X[amb])
    alpha = fit_alpha(X[amb & (tot > 0)], p)
    knee, infl = knee_point(tot, lower)
    test = np.where((tot > lower) & (tot < infl))[0]
    obs = dm_loglik(X[test], p, alpha)
    tmax = int(tot[test].max()) if len(test) else 0
    sims = simulate(p, alpha, tmax, np.unique(tot[test]))
    pv = np.array([(1 + np.searchsorted(sims[int(T)], o, side="right")) / (N_SIM + 1)
                   for T, o in zip(tot[test], obs)])
    # BH
    n = len(pv)
    order = np.argsort(pv)
    q = np.empty(n)
    q[order] = np.minimum.accumulate((pv[order] * n / np.arange(1, n + 1))[::-1])[::-1]
    is_cell = tot >= infl  # retain = inflection (see module docstring)
    is_cell[test[q <= FDR]] = True
    return is_cell, dict(alpha=alpha, knee=knee, inflection=infl, n_tested=n, n_sig=int((q <= FDR).sum()),
                         n_knee=int((tot >= knee).sum()), n_infl=int((tot >= infl).sum()))


def read_raw(d):
    X = scipy.io.mmread(str(d / "matrix.mtx.gz")).T.tocsr().astype(np.float32)
    feat = pd.read_csv(d / "features.tsv.gz", sep="\t", header=None)
    bc = pd.read_csv(d / "barcodes.tsv.gz", sep="\t", header=None)[0].values
    return X, feat, bc


def supp():
    x = pd.read_excel(SUPP, sheet_name=None, header=None)
    si = x["Sample information"].iloc[1:].copy()
    si.columns = si.iloc[0]
    si = si.iloc[1:].set_index("sample_id")
    sn = x["Single-nucleus RNA-seq"].iloc[1:].copy()
    sn.columns = sn.iloc[0]
    sn = sn.iloc[1:].set_index("sample_id")
    return si, sn


def load():
    si, sn = supp()
    parts = []
    for d in sorted(D.glob("GSM*/*_raw_feature_bc_matrix")):
        gsm, sid = re.match(r"(GSM\d+)_(\w+)", d.parent.name).groups()
        donor = re.sub(r"^C(\d\d)$", r"C0\1", sid)
        X, feat, bc = read_raw(d)
        is_cell, info = empty_drops(X)
        g = (feat[2] == "Gene Expression").values
        a = ad.AnnData(X[is_cell][:, g], obs=pd.DataFrame(index=bc[is_cell]),
                       var=pd.DataFrame({"gene_ids": feat[0].values[g]}, index=feat[1].values[g]))
        a.var_names_make_unique()
        Xc = a.X.tocsr()
        ng = np.diff(Xc.indptr)
        tt = np.asarray(Xc.sum(1)).ravel()
        mt = np.asarray(Xc[:, a.var_names.str.startswith("MT-")].sum(1)).ravel()
        mtp = 100 * mt / np.maximum(tt, 1)
        keep = (ng >= 200) & (mtp < 10)
        n5 = int(((ng >= 200) & (mtp < 5)).sum())
        a = a[keep].copy()
        paper = int(sn.loc[donor, "Total.Viable.Nuclei"])
        print(f"{donor}: emptyDrops cells {int(is_cell.sum())} (knee {info['knee']:.0f} UMI -> {info['n_knee']}, "
              f"inflection {info['inflection']:.0f} UMI -> {info['n_infl']} retained, "
              f"tested {info['n_tested']}, sig {info['n_sig']}, alpha {info['alpha']:.0f}); after QC {a.n_obs} "
              f"(with <5% mito {n5}); paper {paper}", flush=True)
        m = si.loc[donor]
        a.obs_names = donor + "_" + a.obs_names
        a.obs["sample"] = donor
        a.obs["gsm"] = gsm
        a.obs["donor"] = donor
        a.obs["region"] = "cingulate"
        grp = m["disease_group"]
        a.obs["disease"] = {"Control": "control", "PD": "PD", "PDD": "PDD", "DLB": "DLB"}[grp]
        a.obs["disease_raw"] = grp
        a.obs["sex"] = m["Sex"]
        a.obs["age"] = float(m["AoD"])
        a.obs["pmi"] = pd.to_numeric(m["PMI"], errors="coerce")
        a.obs["braak_lb"] = pd.to_numeric(m["aSN"], errors="coerce")
        a.obs["braak_tau"] = pd.to_numeric(m["TAU"], errors="coerce")
        a.obs["thal"] = pd.to_numeric(m["thal AB"], errors="coerce")
        a.obs["disease_duration"] = pd.to_numeric(m["DD"], errors="coerce")
        a.obs["acg_asn_score"] = pd.to_numeric(m["aCG aSN score"], errors="coerce")
        a.obs["lpc_category"] = str(m["lpc_category"])
        a.obs["author_n_nuclei"] = paper
        a.obs["author_cell_type"] = "unknown"
        parts.append(a)
    return ad.concat(parts, join="inner", merge="same")

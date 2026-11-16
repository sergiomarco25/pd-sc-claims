"""Meta-analysis helpers shared by the biology syntheses (composition, gene-level DE, programmes, stage, genetics)."""
import numpy as np
from scipy import stats


def fixed(est, se):
    """Inverse-variance fixed-effect combination (used within a dataset across its regions)."""
    est, se = np.asarray(est, float), np.asarray(se, float)
    ok = np.isfinite(est) & np.isfinite(se) & (se > 0)
    if not ok.any():
        return np.nan, np.nan
    w = 1 / se[ok] ** 2
    return float((w * est[ok]).sum() / w.sum()), float(np.sqrt(1 / w.sum()))


def dersimonian_laird(est, se):
    """Random-effects meta-analysis (DerSimonian-Laird). Returns dict: est, se, z, p, tau2, I2, Q, k."""
    est, se = np.asarray(est, float), np.asarray(se, float)
    ok = np.isfinite(est) & np.isfinite(se) & (se > 0)
    est, se = est[ok], se[ok]
    k = len(est)
    if k == 0:
        return {"est": np.nan, "se": np.nan, "z": np.nan, "p": np.nan, "tau2": np.nan, "I2": np.nan, "Q": np.nan, "k": 0}
    w = 1 / se ** 2
    mu_f = (w * est).sum() / w.sum()
    q = float((w * (est - mu_f) ** 2).sum())
    c = w.sum() - (w ** 2).sum() / w.sum()
    tau2 = max(0.0, (q - (k - 1)) / c) if k > 1 and c > 0 else 0.0
    wr = 1 / (se ** 2 + tau2)
    mu = float((wr * est).sum() / wr.sum())
    s = float(np.sqrt(1 / wr.sum()))
    z = mu / s
    i2 = max(0.0, (q - (k - 1)) / q) if k > 1 and q > 0 else 0.0
    return {"est": mu, "se": s, "z": z, "p": float(2 * stats.norm.sf(abs(z))), "tau2": tau2, "I2": i2, "Q": q, "k": k}


def bh(p):
    """Benjamini-Hochberg adjusted p-values (NaN kept)."""
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = np.isfinite(p)
    if ok.sum():
        pv = p[ok]
        o = np.argsort(pv)
        r = pv[o] * len(pv) / (np.arange(len(pv)) + 1)
        r = np.minimum.accumulate(r[::-1])[::-1].clip(0, 1)
        q = np.empty_like(pv)
        q[o] = r
        out[ok] = q
    return out

"""Low-memory readers for large h5ad files (shared by the build loaders).

`read_csr` streams a CSR matrix (X, layers/<name> or raw/X) from an h5ad file row-chunk by row-chunk, keeping only the
requested rows and columns, dropping explicit zeros, and storing float32 values with int64 indices. This avoids loading
float64 matrices and unwanted cells into memory.
"""
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

try:
    from anndata.io import read_elem
except ImportError:  # older anndata
    from anndata.experimental import read_elem


def read_df(path, key="obs"):
    """Read obs / var (or raw/var) as a DataFrame without touching the matrix."""
    with h5py.File(path, "r") as h:
        return read_elem(h[key])


def read_csr(path, key="X", rows=None, cols=None, chunk=50_000, check_int=True):
    """Stream a CSR matrix from an h5ad file.

    rows: boolean mask or sorted integer index of rows to keep (None = all).
    cols: integer index of columns to keep, in the desired output order (None = all).
    """
    with h5py.File(path, "r") as h:
        g = h[key]
        enc = g.attrs.get("encoding-type", "")
        if isinstance(enc, bytes):
            enc = enc.decode()
        assert enc == "csr_matrix", f"{path}:{key} is {enc}, expected csr_matrix"
        n_rows, n_cols = g.attrs["shape"]
        indptr = g["indptr"][:]
        if rows is None:
            keep = np.ones(n_rows, bool)
        else:
            rows = np.asarray(rows)
            keep = rows if rows.dtype == bool else np.isin(np.arange(n_rows), rows)
        out_cols = n_cols if cols is None else len(cols)
        blocks = []
        for s in range(0, n_rows, chunk):
            e = min(s + chunk, n_rows)
            k = keep[s:e]
            if not k.any():
                continue
            a, b = indptr[s], indptr[e]
            data = g["data"][a:b]
            ind = g["indices"][a:b]
            m = sp.csr_matrix((data, ind, indptr[s:e + 1] - a), shape=(e - s, n_cols))
            m = m[np.where(k)[0]]
            if cols is not None:
                m = m[:, cols]
            m.eliminate_zeros()
            if check_int and m.nnz:
                probe = m.data[: 100_000]
                assert np.allclose(probe, np.round(probe)), f"{path}:{key} rows {s}-{e} are not integer counts"
            m = sp.csr_matrix(m, dtype=np.float32)
            blocks.append(m)
    if not blocks:
        return sp.csr_matrix((0, out_cols), dtype=np.float32)
    return _vstack_csr(blocks, out_cols)


def _vstack_csr(blocks, n_cols):
    """vstack CSR blocks with int64 indptr (scipy may pick int32 and overflow for >2^31 nnz)."""
    nnz = int(sum(int(b.nnz) for b in blocks))
    n = int(sum(b.shape[0] for b in blocks))
    idt = np.int64 if nnz > 2**31 - 1 else np.int32
    data = np.empty(nnz, np.float32)
    indices = np.empty(nnz, idt)
    indptr = np.empty(n + 1, idt)
    indptr[0] = 0
    p = r = 0
    blocks.reverse()
    while blocks:  # pop as we go so the blocks are freed while the combined matrix fills
        b = blocks.pop()
        k = int(b.nnz)
        data[p:p + k] = b.data
        indices[p:p + k] = b.indices
        indptr[r + 1:r + 1 + b.shape[0]] = b.indptr[1:].astype(np.int64) + p
        p += k
        r += b.shape[0]
    return sp.csr_matrix((data, indices, indptr), shape=(n, n_cols))


def to_float(x):
    """'>=90' / '89+' / '90+' -> 90; unparsable -> NaN."""
    if pd.isna(x):
        return np.nan
    s = str(x).strip().replace(">=", "").replace("+", "").replace(">", "")
    try:
        return float(s)
    except ValueError:
        return np.nan

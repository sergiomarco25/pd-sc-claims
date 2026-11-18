"""Low-memory readers for large h5ad files, for the field's dataset loaders (claimstore engine).

`read_csr` streams a CSR matrix (X, layers/<name> or raw/X) from an h5ad file row-chunk by row-chunk, keeping only the
requested rows and columns, dropping explicit zeros, and storing float32 values (int64 indices above 2^31 values).
This avoids loading float64 matrices and unwanted cells into memory. `read_csr_many` does the same for several files
(or several matrices of one file) into one matrix, filled in place.
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
    return read_csr_many([(path, key, rows, cols)], chunk=chunk, check_int=check_int)


def read_csr_many(sources, chunk=50_000, check_int=True):
    """Stream the rows of one or more CSR matrices into a single matrix, in source order.

    sources: (path, key[, rows[, cols]]) per matrix, rows and cols as in `read_csr`; every source must give the same
    number of output columns (cols puts a source's columns in the output order). The result is filled in place, so
    memory stays at the final matrix (12 bytes per stored value above 2^31 values, else 8) instead of the chunks plus a
    stacked copy; several files (regions or batches of one dataset) are read without stacking matrices afterwards.
    """
    plan, n_out, nnz_max, n = [], None, 0, 0
    for src in sources:
        path, key, rows, cols = (tuple(src) + (None, None))[:4]
        with h5py.File(path, "r") as h:
            g = h[key]
            enc = g.attrs.get("encoding-type", "")
            if isinstance(enc, bytes):
                enc = enc.decode()
            assert enc == "csr_matrix", f"{path}:{key} is {enc}, expected csr_matrix"
            n_rows, n_cols = (int(x) for x in g.attrs["shape"])
            indptr = g["indptr"][:]
        if rows is None:
            keep = np.ones(n_rows, bool)
        else:
            rows = np.asarray(rows)
            keep = rows if rows.dtype == bool else np.isin(np.arange(n_rows), rows)
        oc = n_cols if cols is None else len(cols)
        n_out = oc if n_out is None else n_out
        assert oc == n_out, f"{path}:{key} gives {oc} columns, the first source {n_out}"
        nnz_max += int(np.diff(indptr)[keep].sum())  # upper bound: before dropping columns and explicit zeros
        n += int(keep.sum())
        plan.append((path, key, keep, cols, indptr, n_cols))
    if n == 0:
        return sp.csr_matrix((0, n_out or 0), dtype=np.float32)
    idt = np.int64 if nnz_max > 2**31 - 1 else np.int32
    data = np.empty(nnz_max, np.float32)
    indices = np.empty(nnz_max, idt)
    out_ptr = np.zeros(n + 1, idt)
    p = r = 0
    for path, key, keep, cols, indptr, n_cols in plan:
        with h5py.File(path, "r") as h:
            g = h[key]
            for s in range(0, len(keep), chunk):
                e = min(s + chunk, len(keep))
                k = keep[s:e]
                if not k.any():
                    continue
                a, b = indptr[s], indptr[e]
                m = sp.csr_matrix((g["data"][a:b], g["indices"][a:b], indptr[s:e + 1] - a), shape=(e - s, n_cols))
                m = m[np.where(k)[0]]
                if cols is not None:
                    m = m[:, cols]
                m.eliminate_zeros()
                if check_int and m.nnz:
                    probe = m.data[: 100_000]
                    assert np.allclose(probe, np.round(probe)), f"{path}:{key} rows {s}-{e} are not integer counts"
                q = int(m.nnz)
                data[p:p + q] = m.data
                indices[p:p + q] = m.indices
                out_ptr[r + 1:r + 1 + m.shape[0]] = m.indptr[1:].astype(np.int64) + p
                p += q
                r += m.shape[0]
    if p < 0.75 * nnz_max:  # many values dropped (column subset): release the unused part
        data, indices = data[:p].copy(), indices[:p].copy()
    else:
        data, indices = data[:p], indices[:p]
    return sp.csr_matrix((data, indices, out_ptr), shape=(n, n_out))


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

"""The engine's low-memory h5ad readers (claimstore.h5io) return what scipy slicing of the full matrix returns."""
import sys
from pathlib import Path

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("anndata")
sp = pytest.importorskip("scipy.sparse")

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from claimstore import h5io  # noqa: E402


def counts(n_rows, n_cols, seed):
    rng = np.random.default_rng(seed)
    m = sp.random(n_rows, n_cols, density=0.2, format="csr", random_state=seed, dtype=np.float32)
    m.data = rng.integers(1, 50, m.nnz).astype(np.float32)
    m.data[::17] = 0  # explicit zeros, as some files carry
    return m


def write_csr(path, m, key="X"):
    with h5py.File(path, "a") as h:
        g = h.create_group(key)
        g.attrs["encoding-type"], g.attrs["encoding-version"] = "csr_matrix", "0.1.0"
        g.attrs["shape"] = np.array(m.shape)
        g["data"], g["indices"], g["indptr"] = m.data, m.indices.astype(np.int32), m.indptr.astype(np.int64)


def same(a, b):
    b = sp.csr_matrix(b, copy=True)
    b.eliminate_zeros()
    assert a.shape == b.shape and a.dtype == np.float32
    assert (a != b).nnz == 0
    assert a.nnz == b.nnz  # no explicit zeros left


def test_read_csr_rows_cols_chunks(tmp_path):
    m = counts(300, 40, 0)
    f = tmp_path / "a.h5ad"
    write_csr(f, m)
    write_csr(f, m * 2, key="raw/X")
    same(h5io.read_csr(f), m)
    same(h5io.read_csr(f, "raw/X", chunk=7), m * 2)
    rows = np.zeros(300, bool)
    rows[::3] = True
    cols = np.array([5, 0, 39, 12, 12])
    same(h5io.read_csr(f, rows=rows, cols=cols, chunk=11), m[np.where(rows)[0]][:, cols])
    same(h5io.read_csr(f, rows=np.array([3, 4, 250]), chunk=50), m[[3, 4, 250]])
    assert h5io.read_csr(f, rows=np.zeros(300, bool)).shape == (0, 40)


def test_read_csr_many_stacks_sources_in_order(tmp_path):
    a, b = counts(120, 30, 1), counts(80, 25, 2)
    fa, fb = tmp_path / "a.h5ad", tmp_path / "b.h5ad"
    write_csr(fa, a)
    write_csr(fb, b)
    rows = np.arange(120) % 4 != 0
    ca, cb = np.arange(20), np.arange(5, 25)
    got = h5io.read_csr_many([(fa, "X", rows, ca), (fb, "X", None, cb)], chunk=9)
    same(got, sp.vstack([a[np.where(rows)[0]][:, ca], b[:, cb]]))
    same(h5io.read_csr_many([(fa, "X")]), a)
    with pytest.raises(AssertionError, match="columns"):
        h5io.read_csr_many([(fa, "X"), (fb, "X")])


def test_read_csr_rejects_non_counts(tmp_path):
    m = counts(50, 10, 3)
    m.data = m.data + 0.5
    f = tmp_path / "c.h5ad"
    write_csr(f, m)
    with pytest.raises(AssertionError, match="integer counts"):
        h5io.read_csr(f)
    assert h5io.read_csr(f, check_int=False).nnz == m.nnz

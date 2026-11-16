"""Check every downloaded file and write inventory/download_manifest.tsv (dataset, file, bytes, md5, check, status).

Usage: python -I scripts/04_verify_downloads.py
Checks: .tar → member list readable; .gz → full decompression; .h5/.h5ad → HDF5 opens and top-level keys listed;
.rds/.xlsx/.json → non-empty. Each file is also compared against the size the source server reports, when it reports one.
"""
from pathlib import Path
import gzip
import hashlib
import sys
import tarfile
import urllib.request
import h5py
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))  # noqa: E402
from pdrep import config as cfg  # noqa: E402

RAW = cfg.RAW
OUT = Path(__file__).resolve().parents[1] / "inventory" / "download_manifest.tsv"
GEO = "https://ftp.ncbi.nlm.nih.gov/geo"


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 24), b""):
            h.update(chunk)
    return h.hexdigest()


def remote_size(name):
    """Content-Length of the GEO file this local file came from (None if not a GEO file or no answer)."""
    acc = name.split("_")[0]
    if acc.startswith("GSE"):
        sub = "soft" if name.endswith("_family.soft.gz") else "suppl"
        url = f"{GEO}/series/{acc[:-3]}nnn/{acc}/{sub}/{name}"
    elif acc.startswith("GSM"):
        url = f"{GEO}/samples/{acc[:-3]}nnn/{acc}/suppl/{name}"
    else:
        return None
    try:
        with urllib.request.urlopen(urllib.request.Request(url, method="HEAD"), timeout=60) as r:
            return int(r.headers["Content-Length"])
    except Exception:
        return None


def check(path):
    n = path.name
    if n.endswith(".tar"):
        with tarfile.open(path) as t:
            return f"tar: {len(t.getmembers())} members"
    if n.endswith(".gz") and not n.endswith(".tar.gz"):
        size = 0
        with gzip.open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 24), b""):
                size += len(chunk)
        return f"gzip ok: {size / 1e6:.1f} MB uncompressed"
    if n.endswith((".h5", ".h5ad")):
        with h5py.File(path, "r") as f:
            return "hdf5 keys: " + ",".join(list(f.keys())[:6])
    return "non-empty" if path.stat().st_size > 0 else "EMPTY"


def main():
    rows = []
    for path in sorted(p for p in RAW.rglob("*") if p.is_file()):
        size = path.stat().st_size
        try:
            msg, status = check(path), "ok"
        except Exception as e:  # report and continue
            msg, status = f"{type(e).__name__}: {e}", "FAIL"
        rs = remote_size(path.name)
        if rs is not None and rs != size:
            status, msg = "FAIL", msg + f"; size {size} != remote {rs}"
        rows.append({"dataset": path.relative_to(RAW).parts[0], "file": str(path.relative_to(RAW)), "bytes": size,
                     "remote_bytes": rs if rs is not None else "", "md5": md5(path), "check": msg, "status": status})
        print(status, rows[-1]["file"], msg, flush=True)
    df = pd.DataFrame(rows)
    df.to_csv(OUT, sep="\t", index=False)
    print(df.groupby("dataset").agg(files=("file", "size"), GB=("bytes", lambda b: round(b.sum() / 1e9, 2)),
                                    failed=("status", lambda s: int((s != "ok").sum()))).to_string())
    sys.exit(int((df.status != "ok").any()))


if __name__ == "__main__":
    main()

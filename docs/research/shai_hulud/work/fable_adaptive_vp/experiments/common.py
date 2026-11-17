"""Shared helpers for the experiment drivers: import paths, metadata, resumable
JSON output, resource accounting.

Every result file carries: command line, package versions, platform, cpu
count, thread settings, a hash of src/shadapt/*.py and of the vendored
baseline, and a time stamp.  Drivers save after every arm and skip arms that
are already present (resume after interruption)."""
from __future__ import annotations

import hashlib
import json
import os
import platform
import resource
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SRC = os.path.join(ROOT, "src")
VENDOR_SRC = os.path.join(ROOT, "vendor", "fable_breakthrough", "src")
for _p in (SRC, VENDOR_SRC):
    if _p not in sys.path:
        sys.path.insert(0, _p)
RESULTS = os.path.join(ROOT, "results")
# large binary controls (reference snapshots) are regenerated, not committed
DATA = os.environ.get("SHADAPT_DATA", os.path.join(ROOT, "results", "data"))


def _hash_dir(d):
    h = hashlib.sha256()
    for f in sorted(os.listdir(d)):
        if f.endswith(".py"):
            h.update(f.encode())
            h.update(open(os.path.join(d, f), "rb").read())
    return h.hexdigest()[:16]


def source_hash():
    return _hash_dir(os.path.join(SRC, "shadapt"))


def vendor_hash():
    return _hash_dir(os.path.join(VENDOR_SRC, "shsheet"))


def metadata(script, **extra):
    import numpy
    import scipy
    try:
        import numba
        nb = numba.__version__
        nthreads = numba.get_num_threads()
    except Exception:  # pragma: no cover
        nb, nthreads = None, None
    m = dict(script=os.path.basename(script), command="python3 " + " ".join(sys.argv), python=sys.version.split()[0],
             numpy=numpy.__version__, scipy=scipy.__version__, numba=nb, numba_threads=nthreads,
             platform=platform.platform(), cpu_count=os.cpu_count(),
             env_threads={k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMBA_NUM_THREADS")},
             source_sha256_16=source_hash(), vendor_sha256_16=vendor_hash(), timestamp=time.strftime("%Y-%m-%dT%H:%M:%S"))
    m.update(extra)
    return m


def save(name, obj):
    os.makedirs(RESULTS, exist_ok=True)
    tmp = os.path.join(RESULTS, name + ".tmp")
    with open(tmp, "w") as f:
        json.dump(obj, f, indent=1, default=float)
    os.replace(tmp, os.path.join(RESULTS, name))


def load_partial(name, fresh):
    """Resume support: previously saved dict if present (meta history kept)."""
    path = os.path.join(RESULTS, name)
    if os.path.exists(path):
        with open(path) as f:
            old = json.load(f)
        old.setdefault("meta_history", []).append(old.get("meta"))
        old["meta"] = fresh["meta"]
        for k, v in fresh.items():
            old.setdefault(k, v)
        return old
    return fresh


def peak_rss_mb():
    """Peak resident set size of this process in MiB (ru_maxrss is KiB on Linux)."""
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024.0


def current_rss_mb():
    with open("/proc/self/statm") as f:
        pages = int(f.read().split()[1])
    return pages * os.sysconf("SC_PAGE_SIZE") / 2 ** 20

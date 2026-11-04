"""Thin Python wrapper around the C++ BP+OSD core (build it with `python cpp/build.py`).

    dec = FastBpOsd(HZ, max_iter=100, ms_scaling=0.625, osd_order=0)      # -1 BP only, 0 OSD-0, >0 OSD-CS(order)
    corr, converged, syndrome_ok, iters = dec.decode_batch(syndromes, priors=p_or_P, threads=8)

`priors` are error probabilities, either shape (n,) shared or (B, n) per shot (GKP soft information)."""
import numpy as np
import scipy.sparse as sp

try:
    from . import _fastbp
except ImportError:                       # extension not built
    _fastbp = None

__all__ = ["FastBpOsd", "fast_available", "probs_to_llr"]


def fast_available():
    return _fastbp is not None


def probs_to_llr(p):
    p = np.clip(np.asarray(p, dtype=np.float64), 1e-12, 1 - 1e-12)
    return np.log((1 - p) / p)


class FastBpOsd:
    def __init__(self, H, max_iter=100, ms_scaling=0.625, osd_order=-1):
        if _fastbp is None:
            raise ImportError("C++ extension not built: run `python cpp/build.py` (needs a C++17 compiler and pybind11)")
        Hc = sp.csr_matrix(H, dtype=np.int32); Hc.sort_indices()
        self.m, self.n = Hc.shape
        self._core = _fastbp.FastBpOsd(self.m, self.n, Hc.indptr.astype(np.int32), Hc.indices.astype(np.int32), max_iter, ms_scaling, osd_order)

    def decode_batch(self, syndromes, priors, threads=1):
        S = np.ascontiguousarray(np.atleast_2d(syndromes), dtype=np.uint8)
        return self._core.decode_batch(S, probs_to_llr(priors), int(threads))

    @property
    def has_openmp(self):
        return bool(_fastbp.openmp)

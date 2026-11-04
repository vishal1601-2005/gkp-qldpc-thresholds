"""Export a decoding workload (code-capacity, per-shot GKP soft priors or hard priors) for cpp/bench/bench_native.
Usage: python scripts/export_instance.py n sigma shots mode(hard|soft) osd_order out.bin"""
import sys, numpy as np, scipy.sparse as sp
from gkp_qldpc import make_good_code, sample_shift_residual, posterior_flip_prob, p_flip_exact, probs_to_llr
n, sigma, B, mode, osd, out = int(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3]), sys.argv[4], int(sys.argv[5]), sys.argv[6]
HX, HZ, k, d = make_good_code(n, 200); H = sp.csr_matrix(HZ, dtype=np.int32); H.sort_indices(); m, N = H.shape
rng = np.random.default_rng(1); r, e = sample_shift_residual(sigma, (B, N), rng)
S = ((e @ HZ.toarray().astype(np.uint8).T) % 2).astype(np.uint8)
P = posterior_flip_prob(r, sigma) if mode == "soft" else np.full((B, N), p_flip_exact(sigma))
with open(out, "wb") as f:
    np.array([m, N, H.nnz, B, 100, osd], dtype=np.int32).tofile(f); np.array([0.625]).tofile(f)
    H.indptr.astype(np.int32).tofile(f); H.indices.astype(np.int32).tofile(f); S.tofile(f); probs_to_llr(P).astype(np.float64).tofile(f)
print(f"wrote {out}: code N={N}, m={m}, d1={d}, {B} shots, {mode} priors, sigma={sigma}")

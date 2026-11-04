"""Apples-to-apples decoder benchmark: ldpc (BpOsdDecoder, per-shot prior updates) vs the C++ core (batch API, OpenMP).
Identical shots, identical priors, identical decoder settings (min-sum, scaling 0.625, 100 iterations, OSD-CS order 7).
Reports per (size, decoder-mode): wall time per shot, speedup, block-failure agreement (z-score), BP-convergence rate.
Note: the ldpc loop includes Python call overhead and per-shot prior updates (required for soft information); for the
pure-C++ cost without any Python see cpp/bench/bench_native.cpp.
Usage: python scripts/bench_decoders.py [sizes=16,24,32] [shots=600] [threads=1,4,8] [out=bench_results.csv]
(Build the extension first:  python cpp/build.py)"""
import sys, time, csv, numpy as np
from ldpc import BpOsdDecoder
from gkp_qldpc import (make_good_code, kernel_dense, logical_fail_batch, sample_shift_residual, posterior_flip_prob,
                       p_flip_exact, FastBpOsd, fast_available)

sizes = [int(x) for x in sys.argv[1].split(",")] if len(sys.argv) > 1 else [16, 24, 32]
shots = int(sys.argv[2]) if len(sys.argv) > 2 else 600
threads = [int(x) for x in sys.argv[3].split(",")] if len(sys.argv) > 3 else [1]
out = sys.argv[4] if len(sys.argv) > 4 else "bench_results.csv"
assert fast_available(), "run `python cpp/build.py` first"
POINTS = (("hard", 0.47), ("hard", 0.50), ("soft", 0.52), ("soft", 0.56))
rows = []
print(f"{'N':>5} {'mode':>5} {'sigma':>6} {'ldpc ms/shot':>13} " + " ".join(f"{'C++ x'+str(t)+' ms':>11} {'speedup':>8}" for t in threads) + f" {'fail ldpc':>9} {'fail C++':>9} {'z':>6} {'BP conv':>8}")
for n in sizes:
    HX, HZ, k, d = make_good_code(n, 200); N = HZ.shape[1]; Kd = kernel_dense(HX); HZd = HZ.toarray().astype(np.uint8)
    fast = FastBpOsd(HZ, 100, 0.625, 7)
    for mode, sg in POINTS:
        rng = np.random.default_rng(n * 1000 + int(sg * 100))
        r, e = sample_shift_residual(sg, (shots, N), rng); S = ((e @ HZd.T) % 2).astype(np.uint8)
        Q = posterior_flip_prob(r, sg) if mode == "soft" else np.full((shots, N), p_flip_exact(sg))
        ld = BpOsdDecoder(HZ, error_rate=0.05, max_iter=100, bp_method="minimum_sum", ms_scaling_factor=0.625, osd_method="osd_cs", osd_order=7)
        Rl = np.zeros((shots, N), np.uint8); t0 = time.perf_counter()
        for b in range(shots):
            ld.update_channel_probs(Q[b]); Rl[b] = ld.decode(S[b])
        t_ld = (time.perf_counter() - t0) / shots
        times = []
        for T in threads:
            best = 1e9
            for _ in range(2):
                t0 = time.perf_counter(); C, cv, ok, it = fast.decode_batch(S, Q, threads=T); best = min(best, (time.perf_counter() - t0) / shots)
            times.append(best)
        fl, ff = logical_fail_batch(e ^ Rl, Kd), logical_fail_batch(e ^ C, Kd)
        se = np.sqrt(fl.mean() * (1 - fl.mean()) / shots + ff.mean() * (1 - ff.mean()) / shots) or 1.0
        z = (ff.mean() - fl.mean()) / se
        print(f"{N:>5} {mode:>5} {sg:>6.2f} {t_ld*1e3:>13.3f} " + " ".join(f"{t*1e3:>11.3f} {t_ld/t:>7.1f}x" for t in times) + f" {fl.mean():>9.3f} {ff.mean():>9.3f} {z:>+6.2f} {cv.mean():>8.3f}", flush=True)
        rows.append(dict(N=N, mode=mode, sigma=sg, ldpc_s_per_shot=t_ld, **{f"cpp_t{T}_s_per_shot": t for T, t in zip(threads, times)},
                         fail_ldpc=float(fl.mean()), fail_cpp=float(ff.mean()), z=float(z), bp_converged=float(cv.mean()), identical=float((C == Rl).all(axis=1).mean())))
with open(out, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
print("wrote", out)

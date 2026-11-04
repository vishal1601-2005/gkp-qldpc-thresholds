"""Code-capacity hard-vs-soft crossings for the bivariate-bicycle family, using the C++ decoder.
Independent cross-check of published BB results (e.g. Yao et al.: about 0.500 hard / 0.556 analog crossing in sigma).
Usage: python scripts/bb_code_capacity.py [shots=6000] [threads=1] [out=bb_code_capacity_results.csv]
CSV columns match the other sweeps: mode(0 hard/1 soft), n(=N qubits), codeseed(0), sigma, shots, fails, syndrome_mismatch, k, d1(literature).
Resumable: finished (mode, N, sigma) points are skipped."""
import sys, os, time, numpy as np
from gkp_qldpc import kernel_dense, code_params, FastBpOsd, run_point_fast, BB_CODES, bb_code

shots = int(sys.argv[1]) if len(sys.argv) > 1 else 6000
threads = int(sys.argv[2]) if len(sys.argv) > 2 else 1
out = sys.argv[3] if len(sys.argv) > 3 else "bb_code_capacity_results.csv"
GRID = {"hard": np.round(np.arange(0.440, 0.5601, 0.015), 4), "soft": np.round(np.arange(0.500, 0.6201, 0.015), 4)}
done = set()
if os.path.exists(out):
    old = np.loadtxt(out, delimiter=",", ndmin=2); done = {(int(r[0]), int(r[1]), round(r[3], 4)) for r in old}
t0 = time.time()
for name in BB_CODES:
    HX, HZ, d = bb_code(name); N, k = code_params(HX, HZ); Kd = kernel_dense(HX); dec = FastBpOsd(HZ, 100, 0.625, 7)
    for mode, mi in (("hard", 0), ("soft", 1)):
        for s in GRID[mode]:
            if (mi, N, round(float(s), 4)) in done: continue
            f, sf = run_point_fast(HX, HZ, Kd, float(s), shots, mode, seed=N * 131 + int(s * 1e4) + mi, decoder=dec, threads=threads)
            with open(out, "a") as fh: fh.write(f"{mi},{N},0,{s},{shots},{f},{sf},{k},{d}\n")
            print(f"[[{N},{k},{d}]] {mode} sigma={s}: P_block={f/shots:.4f}   ({time.time()-t0:.0f}s)", flush=True)

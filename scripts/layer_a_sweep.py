"""LAYER (a) v2: screened code ensemble (full-rank, distance >= typical), hard vs soft, per-code records.
Fixes the v1 design flaws: (1) the 'code lottery' (a d=2 code at n=16 ruined whole curves), (2) block failure with k growing ~ n^2.
Run:  python layer_a_sweep.py [shots=1500] [ncodes=4] [nproc=auto] [sizes=12,16,20,24]
Output: layer_a_results.csv   columns: mode(0 hard/1 soft), n, codeseed, sigma, shots, fails, syndrome_mismatch, k, d1
Resumable: finished jobs are skipped.  One line is saved per finished job.
"""
import numpy as np, sys, time, os
from multiprocessing import Pool, cpu_count
from gkp_qldpc import *
from gkp_qldpc import *

SIG = {"hard": np.round(np.arange(0.460, 0.5801, 0.015), 4),
       "soft": np.round(np.arange(0.520, 0.6401, 0.015), 4)}
FN = "layer_a_results.csv"
# large sizes only need the crossing region (saves ~40% time); these are subsets of the small-size grids so neighbouring sizes share sigma points
SIG_BIG = {"hard": np.round(np.arange(0.460, 0.5351, 0.015), 4), "soft": np.round(np.arange(0.520, 0.5951, 0.015), 4)}   # v2b: one extra LOW point so the crossing is bracketed (cheap: low sigma = few decoder failures)
def sig_grid(mode, n): return SIG_BIG[mode] if n >= 32 else SIG[mode]

def task(args):
    mode, n, cseed, sigma, shots = args
    HX, HZ, k, d1 = make_good_code(n, cseed); Kd = kernel_dense(HX)
    f, sf = run_point(HX, HZ, Kd, sigma, shots, mode, seed=cseed * 7919 + int(sigma * 1e4) + (1 if mode == "soft" else 0))
    return (0 if mode == "hard" else 1, n, cseed, sigma, shots, f, sf, k, d1)

if __name__ == "__main__":
    shots = int(sys.argv[1]) if len(sys.argv) > 1 else 1500
    ncodes = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else cpu_count()
    sizes = [int(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else [12, 16, 20, 24]
    jobs = [(m, n, 200 + c, s, shots) for m in ("hard", "soft") for n in sizes for c in range(ncodes) for s in sig_grid(m, n)]
    done = set()
    if os.path.exists(FN):
        old = np.loadtxt(FN, delimiter=",", ndmin=2)
        done = {(int(r[0]), int(r[1]), int(r[2]), round(r[3], 4)) for r in old}
    jobs = [j for j in jobs if ((0 if j[0] == "hard" else 1), j[1], j[2], round(j[3], 4)) not in done]
    print(f"{len(jobs)} jobs to run on {nproc} procs", flush=True)
    t0 = time.time()
    with Pool(nproc) as pool:
        for i, r in enumerate(pool.imap_unordered(task, sorted(jobs, key=lambda j: j[1])), 1):
            with open(FN, "a") as f: f.write(",".join(str(x) for x in r) + "\n")
            print(f"{i}/{len(jobs)} done ({'hard' if r[0]==0 else 'soft'} n={r[1]} code={r[2]} sigma={r[3]}), {time.time()-t0:.0f}s", flush=True)

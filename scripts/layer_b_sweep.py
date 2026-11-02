"""LAYER (b) sweep: phenomenological noise (noisy syndrome readout through GKP ancillas), hard vs soft decoding.
Rounds T = dmin_rule(n) (= code distance of the screened codes). ratio = sigma_ancilla / sigma_data.
Run:  python layer_b_sweep.py [shots=1000] [ncodes=4] [nproc=auto] [sizes=12,16,20,24] [ratio=1.0] [shift=0.0]
shift lowers every sigma in the grid by this amount (needed when a larger ancilla ratio pushes the crossing to smaller sigma).
Output: layer_b_r<ratio>_results.csv  (same columns as layer a: mode,n,codeseed,sigma,shots,fails,sync,k,d1). Resumable, saved per job.
"""
import numpy as np, sys, time, os
from multiprocessing import Pool, cpu_count
from gkp_qldpc import *
from gkp_qldpc import *
from gkp_qldpc import *

SIG = {"hard": np.round(np.arange(0.370, 0.4601, 0.015), 4), "soft": np.round(np.arange(0.430, 0.5201, 0.015), 4)}

def task(args):
    mode, n, cseed, sigma, shots, ratio = args
    HX, HZ, k, d1 = make_good_code(n, cseed); Kd = kernel_dense(HX)
    f, sf = run_point_phenom(HX, HZ, Kd, sigma, ratio, dmin_rule(n), shots, mode,
                             seed=cseed * 7919 + int(sigma * 1e4) + (1 if mode == "soft" else 0))
    return (0 if mode == "hard" else 1, n, cseed, sigma, shots, f, sf, k, d1)

if __name__ == "__main__":
    shots = int(sys.argv[1]) if len(sys.argv) > 1 else 1000
    ncodes = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else cpu_count()
    sizes = [int(x) for x in sys.argv[4].split(",")] if len(sys.argv) > 4 else [12, 16, 20, 24]
    ratio = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
    shift = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0
    FN = f"layer_b_r{ratio:g}_results.csv"
    jobs = [(m, n, 200 + c, s, shots, ratio) for m in ("hard", "soft") for n in sizes for c in range(ncodes) for s in np.round(SIG[m] - shift, 4)]
    done = set()
    if os.path.exists(FN):
        old = np.loadtxt(FN, delimiter=",", ndmin=2)
        done = {(int(r[0]), int(r[1]), int(r[2]), round(r[3], 4)) for r in old}
    jobs = [j for j in jobs if ((0 if j[0] == "hard" else 1), j[1], j[2], round(j[3], 4)) not in done]
    print(f"{len(jobs)} jobs to run on {nproc} procs -> {FN}", flush=True)
    t0 = time.time()
    with Pool(nproc) as pool:
        for i, r in enumerate(pool.imap_unordered(task, sorted(jobs, key=lambda j: j[1])), 1):
            with open(FN, "a") as f: f.write(",".join(str(x) for x in r) + "\n")
            print(f"{i}/{len(jobs)} done ({'hard' if r[0]==0 else 'soft'} n={r[1]} code={r[2]} sigma={r[3]}), {time.time()-t0:.0f}s", flush=True)

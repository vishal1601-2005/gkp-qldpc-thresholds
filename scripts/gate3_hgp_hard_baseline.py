"""GATE 3 - decoder-layer baseline: HGP of random (3,4)-regular codes, independent bit flips (hard decision), BP+OSD.
Published reference (NOT a strict equality target, different decoder): Grospellier et al. 2021 report ~7.5% (BP+SSF, ideal syndrome),
described by the authors as 'evidence' / formally a pseudo-threshold.  Our own BP+OSD crossing is reported as a measurement.
PRE-REGISTERED sanity checks (all must hold):
  (i)  below the crossing, failure rate decreases with code size;  (ii) crossing lies in 0.05-0.11 (between BP-SSF and hashing bound 0.110);
  (iii) BP+OSD failure <= plain BP failure (osd_order=0) at the same point.
Run:  python gate3_hgp_hard_baseline.py [shots=2000] [ncodes=3] [nproc=auto]
Outputs: gate3_results.csv
"""
import numpy as np, sys, time, itertools
from multiprocessing import Pool, cpu_count
from gkp_qldpc import *
from gkp_qldpc import *

SIZES = [16, 24, 32]            # classical n; quantum N = n^2 + (3n/4)^2 = 400, 900, 1600
PS = [0.05, 0.06, 0.07, 0.08, 0.09]

def task(args):
    n, cseed, p, shots, osd = args
    HX, HZ = make_code(n, cseed)
    Kd = kernel_dense(HX)
    sigma = sigma_from_p(p)
    f, sf = run_point(HX, HZ, Kd, sigma, shots, "hard", seed=cseed * 1000 + int(p * 1e4), osd_order=osd)
    return (n, cseed, p, shots, f, sf, osd)

if __name__ == "__main__":
    shots = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
    ncodes = int(sys.argv[2]) if len(sys.argv) > 2 else 3
    nproc = int(sys.argv[3]) if len(sys.argv) > 3 else cpu_count()
    jobs = [(n, 100 + c, p, shots, 10) for n in SIZES for c in range(ncodes) for p in PS]
    jobs += [(16, 100, p, shots, 0) for p in PS]         # plain-BP control (osd_order=0)
    t0 = time.time()
    out = []
    with Pool(nproc) as pool:
        for i, r in enumerate(pool.imap_unordered(task, jobs), 1):
            out.append(r); print(f"{i}/{len(jobs)} jobs done, {time.time()-t0:.0f}s", flush=True)
    np.savetxt("gate3_results.csv", np.array(out), delimiter=",", header="n,codeseed,p,shots,fails,syndrome_mismatch,osd_order")
    A = np.array(out)
    print(f"elapsed {time.time()-t0:.0f}s; failure rate (BP+OSD, osd_order=10), averaged over codes:")
    print("   p   " + "  ".join(f"n={n:>3d}(N={n*n+(3*n//4)**2})" for n in SIZES))
    tab = {}
    for p in PS:
        row = []
        for n in SIZES:
            m = (A[:, 0] == n) & (np.isclose(A[:, 2], p)) & (A[:, 6] == 10)
            row.append(A[m, 4].sum() / A[m, 3].sum())
        tab[p] = row; print(f"{p:.3f}  " + "  ".join(f"{x:14.4f}" for x in row))
    m0 = (A[:, 6] == 0); m1 = (A[:, 0] == 16) & (A[:, 6] == 10) & (A[:, 1] == 100)
    print("\nplain BP vs BP+OSD at n=16, code seed 100:")
    for p in PS:
        f0 = A[m0 & np.isclose(A[:, 2], p), 4].sum() / A[m0 & np.isclose(A[:, 2], p), 3].sum()
        f1 = A[m1 & np.isclose(A[:, 2], p), 4].sum() / A[m1 & np.isclose(A[:, 2], p), 3].sum()
        print(f"  p={p:.3f}  BP: {f0:.4f}   BP+OSD: {f1:.4f}   {'ok' if f1 <= f0 + 3*np.sqrt(f0*(1-f0)/shots) else 'CHECK'}")

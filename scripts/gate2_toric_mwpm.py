"""GATE 2 - classical-layer baseline: toric code, independent bit flips, MWPM (pymatching).
Known result: threshold 0.1031 +/- 0.0001 (Wang-Harrington-Preskill 2003).
PRE-REGISTERED PASS: fitted p_c within 0.002 of 0.1031 (absolute) using L = 8,12,16,20 and the quadratic-collapse fit.
Run:  python gate2_toric_mwpm.py [shots_per_point=20000]
"""
import numpy as np, sys, time, pymatching
from gkp_qldpc import toric, kernel_dense, logical_fail_batch, fss_fit

def simulate(L, p, shots, rng, matcher, HZ, Kd, batch=2000):
    N = HZ.shape[1]; fails = 0; done = 0
    HZd = HZ.toarray().astype(np.uint8)
    while done < shots:
        B = min(batch, shots - done)
        E = (rng.random((B, N)) < p).astype(np.uint8)
        S = (E @ HZd.T) % 2
        C = np.array([matcher.decode(s) for s in S], dtype=np.uint8)
        fails += int(logical_fail_batch(E ^ C, Kd).sum()); done += B
    return fails

if __name__ == "__main__":
    shots = int(sys.argv[1]) if len(sys.argv) > 1 else 20000
    Ls = [8, 12, 16, 20]; ps = np.linspace(0.094, 0.112, 7)
    rng = np.random.default_rng(12345)
    rows = []; t0 = time.time()
    for L in Ls:
        HX, HZ = toric(L); Kd = kernel_dense(HX)
        matcher = pymatching.Matching(HZ)
        for p in ps:
            f = simulate(L, p, shots, rng, matcher, HZ, Kd)
            rows.append((L, p, f, shots)); print(f"L={L:2d} p={p:.4f} fail={f/shots:.4f}  [{time.time()-t0:.0f}s]", flush=True)
    A = np.array(rows)
    np.savetxt("gate2_toric_results.csv", A, delimiter=",", header="L,p,fails,shots")
    pc, nu, epc, enu, chi2 = fss_fit(A[:, 0], A[:, 1], A[:, 2], A[:, 3])
    print(f"\nFSS fit: p_c = {pc:.5f} +/- {epc:.5f},  nu = {nu:.2f} +/- {enu:.2f},  reduced chi2 = {chi2:.2f}")
    print("known p_c = 0.1031 +/- 0.0001 ->", "PASS" if abs(pc - 0.1031) < 0.002 else "FAIL", "(tolerance 0.002)")

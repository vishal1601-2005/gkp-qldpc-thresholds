"""Pairwise finite-size crossings from a results CSV, with distance-aware pair classification and loss-space translation.
Usage: python analyze_pairwise.py results.csv   (default data/raw/layer_a_results.csv)
Columns: mode(0 hard/1 soft), n, codeseed, sigma, shots, fails, syndrome_mismatch, k, d1"""
import sys, numpy as np
from gkp_qldpc import pairwise_crossings, s0sq_from_dB, eta_threshold

fn = sys.argv[1] if len(sys.argv) > 1 else "data/raw/layer_a_results.csv"
A = np.loadtxt(fn, delimiter=",", ndmin=2); print("reading", fn)
TAG = {"clean": "CLEAN (one distance per size, distance increases)", "same": "SAME distance (multiplicity-dominated, do NOT quote)",
       "mixed": "MIXED distances inside a size (bracket only)"}
allres = {}
for mode, name in ((0, "hard"), (1, "soft")):
    M = A[A[:, 0] == mode]
    if not len(M): continue
    print(f"\n[{name}]")
    allres[name] = pairwise_crossings(M)
    for r in allres[name]:
        print(f"  n={r['n_a']:2d} vs {r['n_b']:2d}: sigma_x = {r['sigma_x']:.4f} +/- {r['err']:.4f} (crossing in {r['n_found']}/300 bootstraps)  [{TAG[r['category']]}]")
if len(allres) == 2:
    sh = [round(s["sigma_x"] - h["sigma_x"], 3) for h, s in zip(allres["hard"], allres["soft"]) if np.isfinite(h["sigma_x"]) and np.isfinite(s["sigma_x"])]
    print("\nsoft - hard crossing shift per pair:", sh)
for cat in ("clean", "mixed"):
    rows = [(n, r) for n, v in allres.items() for r in v if r["category"] == cat and np.isfinite(r["sigma_x"])]
    if not rows: continue
    print(f"\nthreshold LOSS (%) from {cat.upper()} pairs, 10 dB squeezing: rescale / amp")
    for n, r in rows:
        s0 = s0sq_from_dB(10)
        print(f"  {n:>4} {r['n_a']}v{r['n_b']}  sigma_x={r['sigma_x']:.4f}:  {100*(1-eta_threshold(r['sigma_x'], s0, 'rescale')):5.1f} / {100*(1-eta_threshold(r['sigma_x'], s0, 'amp')):5.1f}")

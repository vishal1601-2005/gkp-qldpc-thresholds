"""GATE 0: threshold loss (%) implied by published outer-code thresholds, under each loss-channel convention.
The 'paper' column reproduces a formula without lattice rescaling (WRONG; kept to show how much it matters)."""
import numpy as np
from gkp_qldpc import s0sq_from_dB, eta_threshold, sigma_from_p

THRESHOLDS = {"hashing 0.1100": 0.110028, "toric/MWPM 0.1031": 0.1031, "HGP BP+SSF 0.075": 0.075, "HGP BP+SSF+synd.noise 0.03": 0.03}
print("threshold LOSS (%) at the flip-probability thresholds (hard-decision), by model")
for dB in (10, 14, 18):
    s0 = s0sq_from_dB(dB)
    for name, p in THRESHOLDS.items():
        sg = sigma_from_p(p)
        row = [100 * (1 - eta_threshold(sg, s0, m)) for m in ("paper", "rescale", "amp")]
        print(f"{dB}dB {name:28s} paper {row[0]:5.1f}  rescale {row[1]:5.1f}  amp {row[2]:5.1f}")

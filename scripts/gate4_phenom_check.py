"""GATE 4: the phenomenological code must reduce to the code-capacity result when T=1 and measurement noise -> 0 (ratio tiny).
PASS: block failure agrees with layer-(a) run_point within 3 binomial sigma, hard AND soft, at two sigma values.
Run: python gate4_phenom_check.py [shots=600]"""
import sys, numpy as np
from gkp_qldpc import *
from gkp_qldpc import *
from gkp_qldpc import *
shots = int(sys.argv[1]) if len(sys.argv) > 1 else 600
HX, HZ, k, d = make_good_code(16, 200); Kd = kernel_dense(HX)
ok = True
for mode, sig in (("hard", 0.47), ("hard", 0.50), ("soft", 0.52), ("soft", 0.55)):
    f1, _ = run_point(HX, HZ, Kd, sig, shots, mode, seed=11, osd_order=7)
    f2, _ = run_point_phenom(HX, HZ, Kd, sig, 1e-3, 1, shots, mode, seed=12, osd_order=7)
    p1, p2 = f1 / shots, f2 / shots; se = np.sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / shots)
    z = (p2 - p1) / max(se, 1e-9); ok &= abs(z) < 3
    print(f"{mode} sigma={sig}: code-capacity P={p1:.3f}  phenom(T=1,no meas noise) P={p2:.3f}  z={z:+.2f}")
print("GATE 4:", "PASS" if ok else "FAIL")

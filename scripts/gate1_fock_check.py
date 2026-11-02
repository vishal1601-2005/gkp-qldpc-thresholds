"""GATE 1: exact Fock-space check of the GKP loss channel against the analytic effective-variance models.
Usage: python gate1_fock_check.py [dB ...]     (default 8 10 12 14 16)
Pass criteria (variance space): 'rescale' model within 1% for every point with p_L>1e-4 and eta<1; the v1 formula must be clearly off.
Rows with p_L < 1e-12 are dominated by Fock truncation / grid precision and are flagged, not judged."""
import sys, numpy as np
from gkp_qldpc import fock, sigma2_eff, p_flip_exact, sigma_from_p

dBs = [float(a) for a in sys.argv[1:]] or [8, 10, 12, 14, 16]
etas = (1.0, 0.99, 0.98, 0.95, 0.9, 0.8)
bad = []; discr = False
print(f"{'dB':>4} {'eta':>5} {'p_Fock':>11} {'p_rescale':>11} {'p_paper':>11} | variance deviation: rescale  paper")
for dB in dBs:
    res = fock.logical_error_vs_loss(dB, etas)
    for eta in etas:
        p, s0 = res[eta]
        vr, vp = sigma2_eff(s0, eta, "rescale"), sigma2_eff(s0, eta, "paper")
        vf = sigma_from_p(p) ** 2 if 1e-15 < p < 0.4 else float("nan")
        dr, dp = (vf - vr) / vr, (vf - vp) / vp
        flag = "  (below 1e-12: numerical floor)" if p < 1e-12 else ""
        print(f"{dB:>4.0f} {eta:>5.2f} {p:>11.3e} {p_flip_exact(np.sqrt(vr)):>11.3e} {p_flip_exact(np.sqrt(vp)):>11.3e} |              {dr:+8.2%} {dp:+8.2%}{flag}")
        if eta < 1 and p > 1e-4 and abs(dr) > 0.01: bad.append((dB, eta, round(dr, 4)))
        if eta <= 0.9 and p > 1e-4 and abs(dp) > 0.05: discr = True
print("\nGATE 1:", "PASS" if (not bad and discr) else "FAIL", "| rescale violations:", bad, "| discriminates vs no-rescale formula:", discr)

"""Code-capacity decoding of X-type errors on an HGP code with BP+OSD (hard vs GKP-analog 'soft' priors)."""
import numpy as np
from ldpc import BpOsdDecoder
from .channel import p_flip_exact, sample_shift_residual, posterior_flip_prob
from .codes import logical_fail_batch

def make_decoder(HZ, p, max_iter=100, osd_order=10, ms_scaling=0.625):
    return BpOsdDecoder(HZ, error_rate=float(p), max_iter=max_iter, bp_method="minimum_sum",
                        ms_scaling_factor=ms_scaling, osd_method="osd_cs", osd_order=osd_order)

def run_point(HX, HZ, Kd, sigma, shots, mode, seed, max_iter=100, osd_order=10):
    rng = np.random.default_rng(seed)
    N = HZ.shape[1]; HZd = HZ.toarray().astype(np.uint8)
    p_hard = p_flip_exact(sigma)
    dec = make_decoder(HZ, p_hard, max_iter, osd_order)
    fails = 0; sync_fail = 0
    for _ in range(shots):
        r, e = sample_shift_residual(sigma, N, rng)
        s = (HZd @ e) % 2
        if mode == "soft":
            dec.update_channel_probs(posterior_flip_prob(r, sigma))
        c = dec.decode(s).astype(np.uint8)
        if ((HZd @ c) % 2 != s).any():
            fails += 1; sync_fail += 1; continue
        if logical_fail_batch((e ^ c)[None, :], Kd)[0]:
            fails += 1
    return fails, sync_fail

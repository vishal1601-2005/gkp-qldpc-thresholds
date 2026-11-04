"""Vectorised Monte-Carlo points using the C++ decoder (batch API). Same physics and failure definition as
gkp_qldpc.decoders.run_point and gkp_qldpc.phenom.run_point_phenom, but all shots are sampled at once and decoded by FastBpOsd.

Failure = decoder output does not reproduce the syndrome/detectors OR the residual data error is a nontrivial logical."""
import numpy as np
import scipy.sparse as sp
from .channel import p_flip_exact, sample_shift_residual, posterior_flip_prob
from .codes import logical_fail_batch
from .fastdec import FastBpOsd
from .phenom import spacetime_matrix

__all__ = ["run_point_fast", "run_point_phenom_fast"]


def run_point_fast(HX, HZ, Kd, sigma, shots, mode, seed, max_iter=100, osd_order=7, threads=1, batch=4000, decoder=None):
    rng = np.random.default_rng(seed)
    HZc = sp.csr_matrix(HZ, dtype=np.uint8); N = HZc.shape[1]
    dec = decoder or FastBpOsd(HZ, max_iter, 0.625, osd_order)
    p_hard = p_flip_exact(sigma); fails = sync = 0; done = 0
    while done < shots:
        B = min(batch, shots - done)
        r, e = sample_shift_residual(sigma, (B, N), rng)
        S = np.asarray((HZc @ e.T.astype(np.uint8)).T % 2, dtype=np.uint8)
        P = posterior_flip_prob(r, sigma) if mode == "soft" else np.full(N, p_hard)
        C, conv, ok, it = dec.decode_batch(S, P, threads=threads)
        bad = ok == 0
        fails += int((bad | logical_fail_batch(e ^ C, Kd)).sum()); sync += int(bad.sum()); done += B
    return fails, sync


def run_point_phenom_fast(HX, HZ, Kd, sigma, ratio, T, shots, mode, seed, max_iter=100, osd_order=7, threads=1, batch=500, decoder=None):
    rng = np.random.default_rng(seed)
    HZc = sp.csr_matrix(HZ, dtype=np.uint8); m, N = HZc.shape
    sm = ratio * sigma; pd, pm = p_flip_exact(sigma), p_flip_exact(sm)
    Hst = spacetime_matrix(HZ, T)
    dec = decoder or FastBpOsd(Hst, max_iter, 0.625, osd_order)
    hard_prior = np.concatenate([np.full(N * T, pd), np.full(m * T, pm)])
    fails = sync = 0; done = 0
    while done < shots:
        B = min(batch, shots - done)
        r, e = sample_shift_residual(sigma, (B, T, N), rng)
        ra, ea = sample_shift_residual(sm, (B, T, m), rng)
        E = np.bitwise_xor.reduce(e, axis=1)                                   # (B, N) accumulated data error
        cum = (np.cumsum(e, axis=1) % 2).astype(np.uint8)                      # (B, T, N)
        s_true = np.asarray((HZc @ cum.reshape(B * T, N).T).T % 2, dtype=np.uint8).reshape(B, T, m)
        s_meas = s_true ^ ea
        s_prev = np.concatenate([np.zeros((B, 1, m), np.uint8), s_meas], axis=1)
        s_all = np.concatenate([s_meas, s_true[:, -1:, :]], axis=1)
        det = (s_all ^ s_prev).reshape(B, (T + 1) * m).astype(np.uint8)
        if mode == "soft":
            P = np.concatenate([posterior_flip_prob(r, sigma).reshape(B, T * N), posterior_flip_prob(ra, sm).reshape(B, T * m)], axis=1)
        else:
            P = hard_prior
        C, conv, ok, it = dec.decode_batch(det, P, threads=threads)
        Cdata = np.bitwise_xor.reduce(C[:, :N * T].reshape(B, T, N), axis=1)
        bad = ok == 0
        fails += int((bad | logical_fail_batch(E ^ Cdata, Kd)).sum()); sync += int(bad.sum()); done += B
    return fails, sync

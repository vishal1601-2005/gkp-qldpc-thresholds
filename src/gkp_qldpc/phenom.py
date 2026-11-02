"""LAYER (b): phenomenological noise with GKP analog information.
Model (X-type errors only; Z-type is symmetric for CSS + identical independent noise):
  * T rounds. In each round every data qubit suffers a Gaussian GKP shift of std sigma  -> bit flip + analog residual r (soft info).
  * In each round every check is read out through an ancilla GKP qubit with shift std sigma_m = ratio * sigma -> measurement flip + analog residual.
  * After T noisy rounds, one perfect readout round (standard phenomenological convention).
  * Decoder: BP+OSD on the space-time detector matrix. 'hard': uniform priors p(sigma), p(sigma_m); 'soft': per-variable posteriors from the residuals.
Free modelling parameters (state them in the paper): T (rounds), ratio=sigma_m/sigma, measurement noise independent of check weight.
Failure = detectors not reproduced OR residual data error is a nontrivial logical.
"""
import numpy as np, scipy.sparse as sp
from ldpc import BpOsdDecoder
from .channel import p_flip_exact, sample_shift_residual, posterior_flip_prob
from .codes import logical_fail_batch

def spacetime_matrix(HZ, T):
    """Detectors d_t, t=0..T (T+1 blocks of m rows). Variables: data errors e_1..e_T (N each) then measurement errors m_1..m_T (m each)."""
    HZ = sp.csr_matrix(HZ, dtype=np.int64); m, N = HZ.shape
    data = sp.lil_matrix((m * (T + 1), N * T), dtype=np.int64)
    for t in range(T): data[t * m:(t + 1) * m, t * N:(t + 1) * N] = HZ
    meas = sp.lil_matrix((m * (T + 1), m * T), dtype=np.int64)
    I = sp.identity(m, dtype=np.int64, format="lil")
    for t in range(T):
        meas[t * m:(t + 1) * m, t * m:(t + 1) * m] = I          # m_{t+1} flips detector block t ...
        meas[(t + 1) * m:(t + 2) * m, t * m:(t + 1) * m] = I    # ... and block t+1
    H = sp.hstack([data, meas]).tocsr(); H.data %= 2; H.eliminate_zeros()
    return H.astype(np.uint8)

def run_point_phenom(HX, HZ, Kd, sigma, ratio, T, shots, mode, seed, max_iter=100, osd_order=7):
    rng = np.random.default_rng(seed)
    HZc = sp.csr_matrix(HZ); m, N = HZc.shape
    HZd = HZc.toarray().astype(np.uint8)
    Hst = spacetime_matrix(HZ, T); sm = ratio * sigma
    pd, pm = p_flip_exact(sigma), p_flip_exact(sm)
    prior = np.concatenate([np.full(N * T, pd), np.full(m * T, pm)])
    dec = BpOsdDecoder(Hst, error_rate=float(pd), max_iter=max_iter, bp_method="minimum_sum",
                       ms_scaling_factor=0.625, osd_method="osd_cs", osd_order=osd_order)
    dec.update_channel_probs(prior)
    fails = sync = 0
    for _ in range(shots):
        r, e = sample_shift_residual(sigma, (T, N), rng)             # data shifts per round
        ra, ea = sample_shift_residual(sm, (T, m), rng)              # ancilla shifts per round
        E = np.bitwise_xor.reduce(e, axis=0)                          # accumulated data error
        cum = np.cumsum(e, axis=0) % 2                                # data error after each round
        s_true = (cum @ HZd.T) % 2                                    # (T,m) true syndromes
        s_meas = s_true ^ ea                                          # noisy readouts
        s_prev = np.vstack([np.zeros((1, m), np.uint8), s_meas])      # sigma_0 = 0
        s_all = np.vstack([s_meas, s_true[-1:]])                      # sigma_{T+1} = perfect
        det = (s_all ^ s_prev).astype(np.uint8).reshape(-1)           # (T+1)*m detectors
        if mode == "soft":
            q = np.concatenate([posterior_flip_prob(r, sigma).reshape(-1), posterior_flip_prob(ra, sm).reshape(-1)])
            dec.update_channel_probs(q)
        c = dec.decode(det).astype(np.uint8)
        if ((Hst @ c) % 2 != det).any(): fails += 1; sync += 1; continue
        C = np.bitwise_xor.reduce(c[:N * T].reshape(T, N), axis=0)
        if logical_fail_batch((E ^ C)[None, :], Kd)[0]: fails += 1
    return fails, sync

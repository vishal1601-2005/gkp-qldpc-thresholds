"""C++ decoder validation against the `ldpc` package. These are the checks that were run when the core was written."""
import numpy as np
import pytest
from ldpc import BpDecoder, BpOsdDecoder
from gkp_qldpc import (make_good_code, kernel_dense, logical_fail_batch, sample_shift_residual, posterior_flip_prob, p_flip_exact,
                       FastBpOsd, fast_available, run_point, run_point_fast, run_point_phenom, run_point_phenom_fast)

pytestmark = pytest.mark.skipif(not fast_available(), reason="C++ extension not built (python cpp/build.py)")


@pytest.fixture(scope="module")
def code():
    HX, HZ, k, d = make_good_code(12, 200)
    return HX, HZ, kernel_dense(HX), HZ.toarray().astype(np.uint8)


def _shots(HZd, N, sigma, B, mode, seed):
    rng = np.random.default_rng(seed); r, e = sample_shift_residual(sigma, (B, N), rng)
    S = ((e @ HZd.T) % 2).astype(np.uint8)
    Q = posterior_flip_prob(r, sigma) if mode == "soft" else np.full((B, N), p_flip_exact(sigma))
    return e, S, Q


@pytest.mark.parametrize("sigma,mode", [(0.45, "hard"), (0.52, "hard"), (0.55, "soft")])
def test_bp_stage_is_bit_identical_to_ldpc(code, sigma, mode):
    HX, HZ, Kd, HZd = code; N = HZ.shape[1]; B = 400
    e, S, Q = _shots(HZd, N, sigma, B, mode, 1)
    ref = BpDecoder(HZ, error_rate=0.05, max_iter=100, bp_method="minimum_sum", ms_scaling_factor=0.625)
    R = np.zeros((B, N), np.uint8); conv = np.zeros(B, bool); it = np.zeros(B, int)
    for b in range(B):
        ref.update_channel_probs(Q[b]); R[b] = ref.decode(S[b]); conv[b] = ref.converge; it[b] = ref.iter
    C, cv, ok, iters = FastBpOsd(HZ, 100, 0.625, -1).decode_batch(S, Q)
    assert (C == R).all()
    assert (cv.astype(bool) == conv).all() and (iters == it).all()


def test_osd0_matches_ldpc_and_always_satisfies_syndrome(code):
    HX, HZ, Kd, HZd = code; N = HZ.shape[1]; B = 500
    e, S, Q = _shots(HZd, N, 0.50, B, "hard", 2)
    ld = BpOsdDecoder(HZ, error_rate=0.05, max_iter=100, bp_method="minimum_sum", ms_scaling_factor=0.625, osd_method="osd_0", osd_order=0)
    R = np.zeros((B, N), np.uint8)
    for b in range(B): ld.update_channel_probs(Q[b]); R[b] = ld.decode(S[b])
    C, cv, ok, it = FastBpOsd(HZ, 100, 0.625, 0).decode_batch(S, Q)
    assert ((C @ HZd.T) % 2 == S).all()                       # C++ output always reproduces the syndrome
    assert (C == R).all(axis=1).mean() > 0.99                 # differences are only tie-breaking in the column sort
    fl, ff = logical_fail_batch(e ^ R, Kd).mean(), logical_fail_batch(e ^ C, Kd).mean()
    assert abs(fl - ff) < 0.05


@pytest.mark.slow
@pytest.mark.parametrize("sigma,mode", [(0.50, "hard"), (0.54, "soft")])
def test_osd_cs_failure_rate_matches_ldpc(code, sigma, mode):
    HX, HZ, Kd, HZd = code; N = HZ.shape[1]; B = 1500
    e, S, Q = _shots(HZd, N, sigma, B, mode, 3)
    ld = BpOsdDecoder(HZ, error_rate=0.05, max_iter=100, bp_method="minimum_sum", ms_scaling_factor=0.625, osd_method="osd_cs", osd_order=7)
    R = np.zeros((B, N), np.uint8)
    for b in range(B): ld.update_channel_probs(Q[b]); R[b] = ld.decode(S[b])
    C, cv, ok, it = FastBpOsd(HZ, 100, 0.625, 7).decode_batch(S, Q)
    fl, ff = logical_fail_batch(e ^ R, Kd).mean(), logical_fail_batch(e ^ C, Kd).mean()
    se = np.sqrt(fl * (1 - fl) / B + ff * (1 - ff) / B)
    assert abs(ff - fl) < 3 * se
    assert (C == R).all(axis=1).mean() > 0.85


def test_batch_decoding_is_thread_invariant(code):
    HX, HZ, Kd, HZd = code; N = HZ.shape[1]
    e, S, Q = _shots(HZd, N, 0.5, 200, "soft", 4)
    dec = FastBpOsd(HZ, 100, 0.625, 3)
    a = dec.decode_batch(S, Q, threads=1)[0]; b = dec.decode_batch(S, Q, threads=3)[0]
    assert (a == b).all()


def test_input_validation(code):
    HX, HZ, Kd, HZd = code; dec = FastBpOsd(HZ, 100, 0.625, 0); m, N = HZ.shape
    with pytest.raises(Exception): dec.decode_batch(np.zeros((3, m + 1), np.uint8), np.full(N, 0.05))
    with pytest.raises(Exception): dec.decode_batch(np.zeros((3, m), np.uint8), np.full(N + 1, 0.05))


@pytest.mark.slow
def test_fast_code_capacity_point_matches_python_implementation(code):
    HX, HZ, Kd, HZd = code; n = 2000
    f1, _ = run_point(HX, HZ, Kd, 0.52, n, "soft", seed=1, osd_order=7); f2, _ = run_point_fast(HX, HZ, Kd, 0.52, n, "soft", seed=2)
    p1, p2 = f1 / n, f2 / n; se = np.sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / n)
    assert abs(p2 - p1) < 3 * se


@pytest.mark.slow
def test_fast_phenomenological_point_matches_python_implementation(code):
    HX, HZ, Kd, HZd = code; n = 400
    f1, _ = run_point_phenom(HX, HZ, Kd, 0.45, 1.0, 3, n, "hard", seed=1, osd_order=7)
    f2, _ = run_point_phenom_fast(HX, HZ, Kd, 0.45, 1.0, 3, n, "hard", seed=2)
    p1, p2 = f1 / n, f2 / n; se = np.sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / n) or 1.0
    assert abs(p2 - p1) < 3.5 * se

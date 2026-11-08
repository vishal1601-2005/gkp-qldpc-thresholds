"""Golden-output regression for the C++ decoder: exact outputs (SHA-256 of all corrections) on fixed shots. Hashes were generated with the
original OSD implementation, so any optimisation of the C++ core that changes even one correction fails this test."""
import hashlib
import numpy as np
import pytest
from gkp_qldpc import make_good_code, sample_shift_residual, posterior_flip_prob, p_flip_exact, FastBpOsd, fast_available

pytestmark = pytest.mark.skipif(not fast_available(), reason="C++ extension not built (python cpp/build.py)")
GOLDEN = {('hard', 0): '3dae4acbabe40b34e412534e0e96708503c8535f4b0413f101173f6f4fce67a5', ('hard', 7): 'be66aa2a6acae86e3bce3f0b63f532edc6934baf7107aaf0bc0e7ecb9d2c5c99', ('soft', 0): '99b03e55971e6353870029ef1380329bac978070d2de49527281e9866e05a22b', ('soft', 7): '5374a00d05352993b583376bae4d4d0fed3dc994be6066cc59d0b9bb0fc99238'}


@pytest.mark.parametrize("mode,sigma", [("hard", 0.50), ("soft", 0.55)])
@pytest.mark.parametrize("osd", [0, 7])
def test_outputs_match_golden_hash(mode, sigma, osd):
    HX, HZ, k, d = make_good_code(16, 200); N = HZ.shape[1]; HZd = HZ.toarray().astype(np.uint8)
    rng = np.random.default_rng(123); r, e = sample_shift_residual(sigma, (300, N), rng); S = ((e @ HZd.T) % 2).astype(np.uint8)
    Q = posterior_flip_prob(r, sigma) if mode == "soft" else np.full((300, N), p_flip_exact(sigma))
    C = FastBpOsd(HZ, 100, 0.625, osd).decode_batch(S, Q)[0]
    assert hashlib.sha256(C.tobytes()).hexdigest() == GOLDEN[(mode, osd)]

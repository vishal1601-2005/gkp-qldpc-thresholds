import numpy as np
import pytest
from gkp_qldpc import make_good_code, kernel_dense, run_point, run_point_phenom, spacetime_matrix


@pytest.fixture(scope="module")
def code():
    HX, HZ, k, d = make_good_code(12, 200)
    return HX, HZ, kernel_dense(HX)


def test_soft_information_helps(code):
    HX, HZ, Kd = code
    fh, _ = run_point(HX, HZ, Kd, 0.52, 300, "hard", seed=1)
    fs, _ = run_point(HX, HZ, Kd, 0.52, 300, "soft", seed=1)
    assert fs < fh                      # at sigma=0.52 soft fails far less often (about 0.1 vs 0.6)


def test_failure_increases_with_noise(code):
    HX, HZ, Kd = code
    lo, _ = run_point(HX, HZ, Kd, 0.42, 300, "hard", seed=2)
    hi, _ = run_point(HX, HZ, Kd, 0.52, 300, "hard", seed=2)
    assert hi > lo


def test_spacetime_matrix_shape_and_structure(code):
    HX, HZ, Kd = code
    m, N = HZ.shape; T = 3
    H = spacetime_matrix(HZ, T)
    assert H.shape == (m * (T + 1), N * T + m * T)
    meas = H[:, N * T:].toarray()
    assert (meas.sum(axis=0) == 2).all()        # each measurement fault flips exactly two detectors


@pytest.mark.slow
def test_phenomenological_reduces_to_code_capacity(code):
    """Gate 4: T=1 and (almost) no measurement noise must reproduce code-capacity decoding (soft mode, 3-sigma)."""
    HX, HZ, Kd = code
    n = 600
    f1, _ = run_point(HX, HZ, Kd, 0.55, n, "soft", seed=11, osd_order=7)
    f2, _ = run_point_phenom(HX, HZ, Kd, 0.55, 1e-3, 1, n, "soft", seed=12, osd_order=7)
    p1, p2 = f1 / n, f2 / n
    se = np.sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / n)
    assert abs(p2 - p1) < 3 * se

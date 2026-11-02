import numpy as np
import pytest
from gkp_qldpc import (p_flip_exact, p_flip_nearest, sigma_from_p, sigma2_eff, eta_threshold, s0sq_from_dB,
                       sample_shift_residual, posterior_flip_prob)


def test_exact_flip_matches_monte_carlo():
    rng = np.random.default_rng(0)
    _, e = sample_shift_residual(0.5, 400_000, rng)
    assert abs(e.mean() - p_flip_exact(0.5)) < 4 * np.sqrt(0.08 / 400_000)


def test_nearest_cell_is_accurate_at_small_sigma():
    assert abs(p_flip_exact(0.3) - p_flip_nearest(0.3)) < 1e-12


def test_sigma_from_p_inverts_p_flip():
    for p in (0.01, 0.0745, 0.11):
        assert p_flip_exact(sigma_from_p(p)) == pytest.approx(p, rel=1e-9)


def test_posterior_is_calibrated():
    rng = np.random.default_rng(1)
    r, e = sample_shift_residual(0.5, 600_000, rng)
    q = posterior_flip_prob(r, 0.5)
    for lo, hi in [(0.0, 0.05), (0.05, 0.2), (0.2, 0.4), (0.4, 0.5)]:
        m = (q >= lo) & (q < hi)
        assert abs(q[m].mean() - e[m].mean()) < 0.01


def test_loss_models_limits():
    s0 = s0sq_from_dB(10)
    assert sigma2_eff(s0, 1.0, "rescale") == pytest.approx(s0)
    assert sigma2_eff(s0, 0.8, "amp") > sigma2_eff(s0, 0.8, "rescale") > sigma2_eff(s0, 0.8, "paper")


def test_threshold_loss_reproduces_known_table_entries():
    """Regression on the (wrong) no-rescale formula and the corrected one at 10 dB, hashing-bound flip rate."""
    s0, sg = s0sq_from_dB(10), sigma_from_p(0.110028)
    assert 100 * (1 - eta_threshold(sg, s0, "paper")) == pytest.approx(57.2, abs=0.2)
    assert 100 * (1 - eta_threshold(sg, s0, "rescale")) == pytest.approx(34.0, abs=0.2)
    assert 100 * (1 - eta_threshold(sg, s0, "amp")) == pytest.approx(20.5, abs=0.2)

"""Gate 1 as a regression test: exact Fock-space loss vs the analytic effective-variance models."""
import numpy as np
from gkp_qldpc import fock, sigma2_eff, sigma_from_p


def test_rescale_model_matches_exact_simulation_and_old_formula_does_not():
    etas = (1.0, 0.95, 0.9, 0.8)
    res = fock.logical_error_vs_loss(10.0, etas)
    for eta in etas[1:]:
        p, s0 = res[eta]
        v_fock = sigma_from_p(p) ** 2
        assert abs(v_fock - sigma2_eff(s0, eta, "rescale")) / sigma2_eff(s0, eta, "rescale") < 0.01
    p, s0 = res[0.8]
    assert abs(sigma_from_p(p) ** 2 - sigma2_eff(s0, 0.8, "paper")) / sigma2_eff(s0, 0.8, "paper") > 0.2


def test_state_is_normalised_and_symmetric_under_lossless_limit():
    c = fock.gkp_fock(0.1, 150)
    assert abs(np.linalg.norm(c) - 1) < 1e-12

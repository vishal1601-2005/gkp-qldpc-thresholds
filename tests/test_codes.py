import numpy as np
import pytest
from gkp_qldpc import toric, hgp, rand_regular, code_params, make_good_code, dmin_rule, kernel_dense, logical_fail_batch


def test_toric_code_parameters():
    HX, HZ = toric(6)
    assert code_params(HX, HZ) == (72, 2)
    assert not ((HX @ HZ.T).toarray() % 2).any()


def test_dmin_rule_values():
    assert [dmin_rule(n) for n in (12, 16, 20, 24, 28, 32, 36, 40, 44, 48)] == [4, 4, 6, 6, 8, 8, 10, 10, 12, 12]


def test_screened_code_has_fixed_k_and_required_distance():
    HX, HZ, k, d1 = make_good_code(16, 200)
    assert k == 16 and d1 >= dmin_rule(16)
    assert code_params(HX, HZ)[1] == k
    assert not ((HX @ HZ.T).toarray() % 2).any()


def test_screening_is_deterministic():
    a = make_good_code(12, 200); b = make_good_code(12, 200)
    assert (a[0] != b[0]).nnz == 0 and a[3] == b[3]


def test_logical_check_detects_stabilizer_vs_logical():
    HX, HZ = toric(4)
    Kd = kernel_dense(HX)
    stabilizer = np.asarray(HX.toarray()[0], dtype=np.uint8)[None, :]
    assert not logical_fail_batch(stabilizer, Kd)[0]
    # a non-contractible loop of X errors (one row of horizontal edges in the toric layout) is a logical operator
    L = 4; loop = np.zeros(HX.shape[1], dtype=np.uint8); loop[:L] = 1
    assert logical_fail_batch(loop[None, :], Kd)[0]


def test_bivariate_bicycle_family_parameters():
    from gkp_qldpc.codes import BB_CODES, bb_code
    expected = {"72_12_6": (72, 12), "90_8_10": (90, 8), "108_8_10": (108, 8), "144_12_12": (144, 12), "288_12_18": (288, 12)}
    for name, (N, k) in expected.items():
        HX, HZ, d = bb_code(name)
        assert code_params(HX, HZ) == (N, k)
        assert not ((HX @ HZ.T).toarray() % 2).any()
        assert set(np.asarray(HX.sum(1)).ravel()) == {6} and set(np.asarray(HX.sum(0)).ravel()) == {3}

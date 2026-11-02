import numpy as np
from gkp_qldpc import crossing, classify_pair, pairwise_crossings


def test_crossing_interpolates():
    ca = {0.4: 0.2, 0.5: 0.5, 0.6: 0.8}; cb = {0.4: 0.1, 0.5: 0.5, 0.6: 0.9}
    assert crossing(ca, cb) == np.float64(0.5) or abs(crossing(ca, cb) - 0.5) < 1e-9
    assert np.isnan(crossing({0.4: 0.1, 0.5: 0.2}, {0.4: 0.3, 0.5: 0.4}))


def test_pair_classification():
    assert classify_pair([6], [6]) == "same"
    assert classify_pair([6], [8]) == "clean"
    assert classify_pair([4, 6, 8], [6]) == "mixed"


def test_pairwise_crossings_on_synthetic_data():
    rng = np.random.default_rng(0); rows = []
    for n, d in ((12, 4), (16, 6), (20, 8)):                 # sigmoid curves that cross at 0.5 with sharpness growing in n
        for seed in range(4):
            for s in np.arange(0.40, 0.601, 0.02):
                p = 1 / (1 + np.exp(-(s - 0.5) * n * 2.0)); shots = 4000
                rows.append((0, n, seed, round(s, 4), shots, rng.binomial(shots, p), 0, 4, d))
    out = pairwise_crossings(np.array(rows, float), n_boot=50)
    assert all(abs(r["sigma_x"] - 0.5) < 0.02 for r in out) and all(r["category"] == "clean" for r in out)

# Raw per-job results

Place the raw result files from the production runs here (one line per job; columns
`mode(0 hard/1 soft), n, codeseed, sigma, shots, fails, syndrome_mismatch, k, d1`):

- `layer_a_results.csv`   (code capacity; screened HGP(3,4) codes, n = 12..44)
- `layer_b_r1_results.csv` (phenomenological, ancilla/data noise ratio 1.0, n = 12..24)
- `layer_b_r2_results.csv` (ratio 2.0, n = 12..20)

`python scripts/analyze_pairwise.py data/raw/layer_a_results.csv` regenerates the crossing tables in the README.

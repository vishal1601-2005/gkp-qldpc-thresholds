**Hotspot:** OSD column extraction (building the column bit-vectors of the reduced matrix), which the profile pointed at once BP was ruled out as the easier target (see `cpp/experiments/README.md`).

**Change:** the column bit-vectors are now built with a blocked 64x64 bit-matrix transpose (sequential row reads, sequential writes), instead of one strided bit test per (row, column) pair; see the comments in `cpp/include/gkpbp/decoder.hpp`.

**Result** (N = 1600, hard sigma = 0.50, one thread, five alternating old/new pairs, microseconds per shot):

| | old | new |
|---|---|---|
| range | 4124 - 4165 | 3234 - 3275 |
| mean | 4140 | 3251 |

That is about 21.5 % less time (1.27x), with the same sign in all five pairs and a spread of about 1 %.

**Correctness:** the golden-output hashes in `tests/test_cpp_golden.py` were generated from the original OSD, so any change to even one correction fails the test; the full suite (`pytest -q`) passes with the new code.

**Not adopted:** removing the stored variable-to-check message array gave no gain that could be separated from run-to-run noise (`cpp/experiments/README.md`).

**Untried next candidates** (each must keep the BP stage bit-identical to `ldpc`, so the per-column summation order must stay the same): store messages in variable-major order so the variable-node pass is contiguous, and vectorise across columns of equal degree.

*Note:* absolute timings on this laptop (WSL2) differ from session to session, so only the alternating old/new pair above is comparable; the single-thread line in the Timing section was measured in a different session.

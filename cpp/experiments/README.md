# Experiments (not used by the main build)

`decoder_no_v2c.hpp` drops the stored variable-to-check message array and recomputes `post[var] - c2v[e]` on demand.
The subtraction is the same operation on the same doubles, so the BP stage stays **bit-identical** to `ldpc`
(`tests/test_cpp.py` passes with this header in place). It removes one array (about 24 % of the N = 400 working set) and one pass over the edges.

Measured effect: **none that I could resolve** on a cloud Intel Xeon (48 KiB L1d): N = 900 old 1790-2050 us/shot vs new 1935-2065;
N = 1600 old 4330-4850 vs new 4490-4550 (single core, noisy). It may behave differently on a 32 KiB-L1 core such as Zen 4, where cachegrind
reports a 1.2 % L1 miss rate at N = 400 against 0.1 % with a 48 KiB L1.

To try it: build `cpp/bench/bench_native.cpp` against each header (`-I cpp/include` vs `-I <dir containing gkpbp/decoder.hpp>`),
run `./bench 2>&1 inst24.bin 1 5` several times alternately, and keep it only if the gain is larger than the run-to-run spread.

## Measured on a Ryzen 7 7445HS (32 KiB L1d), WSL2, N = 1600, single thread, 4 alternating runs (author's machine)

| run | original (us/shot) | no-v2c variant |
|---|---|---|
| 1 | 3633.1 | 3591.4 |
| 2 | 3548.6 | 3465.9 |
| 3 | 3526.4 | 3427.5 |
| 4 | 3571.0 | 3438.2 |

Mean 3570 vs 3481 us/shot: about 2.5 % faster, with the same sign in all four pairs but a gain comparable to the run-to-run spread (about 3 %).
Not adopted. The profile pointed at a larger target (OSD column extraction, now optimised in the main header).

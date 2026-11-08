# Profiling the C++ decoder core

Generated from `scripts/profile_all.sh`; raw files are in `docs/profiling/`.

## Machine
```
Linux Vishal 6.18.33.2-microsoft-standard-WSL2 #1 SMP PREEMPT_DYNAMIC Thu Jun 18 21:54:43 UTC 2026 x86_64 x86_64 x86_64 GNU/Linux
CPU(s):                                  12
Model name:                              AMD Ryzen 7 7445HS w/ Radeon 740M Graphics
Thread(s) per core:                      2
L1d cache:                               192 KiB (6 instances)
L2 cache:                                6 MiB (6 instances)
L3 cache:                                16 MiB (1 instance)
g++ (Ubuntu 13.3.0-6ubuntu2~24.04.1) 13.3.0
Python 3.12.3
```

## Timing (single thread, best of 3; N = 900 and N = 1600, hard sigma = 0.50, OSD-CS(7), 300 shots)
```
m=432 n=900 E=3024 B=300 osd_order=7 threads=1 | best of 3: 0.3701 s  1233.6 us/shot  811 shots/s | BP converged 0.017  syndrome ok 1.000  mean iters 98.8
m=768 n=1600 E=5376 B=300 osd_order=7 threads=1 | best of 3: 0.8339 s  2779.5 us/shot  360 shots/s | BP converged 0.023  syndrome ok 1.000  mean iters 98.7
```

## perf
```
46.81%  bench_native_pr  bench_native_prof     [.] main._omp_fn.0
28.56%  bench_native_pr  bench_native_prof     [.] gkpbp::BpOsd::osd_cs(double const*, unsigned char*, int)
22.35%  bench_native_pr  bench_native_prof     [.] gkpbp::BpOsd::osd(unsigned char const*, double const*, unsigned char*, gkpbp::Result&)
0.69%  bench_native_pr  bench_native_prof     [.] void std::__merge_sort_with_buffer<__gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, int*, __gnu_cxx::__ops::_Iter_comp_iter<gkpbp::BpOsd::osd(unsigned char const*, double const*, unsigned char*, gkpbp::Result&)::{lambda(int, int)#1}> >(__gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, __gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, int*, __gnu_cxx::__ops::_Iter_comp_iter<gkpbp::BpOsd::osd(unsigned char const*, double const*, unsigned char*, gkpbp::Result&)::{lambda(int, int)#1}>)
0.54%  bench_native_pr  bench_native_prof     [.] __gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > > std::__move_merge<int*, __gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, __gnu_cxx::__ops::_Iter_comp_iter<gkpbp::BpOsd::osd(unsigned char const*, double const*, unsigned char*, gkpbp::Result&)::{lambda(int, int)#1}> >(int*, int*, int*, int*, __gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, __gnu_cxx::__ops::_Iter_comp_iter<gkpbp::BpOsd::osd(unsigned char const*, double const*, unsigned char*, gkpbp::Result&)::{lambda(int, int)#1}>)
0.32%  bench_native_pr  bench_native_prof     [.] void std::__insertion_sort<__gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, __gnu_cxx::__ops::_Iter_comp_iter<gkpbp::BpOsd::osd(unsigned char const*, double const*, unsigned char*, gkpbp::Result&)::{lambda(int, int)#1}> >(__gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, __gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, __gnu_cxx::__ops::_Iter_comp_iter<gkpbp::BpOsd::osd(unsigned char const*, double const*, unsigned char*, gkpbp::Result&)::{lambda(int, int)#1}>)
0.29%  bench_native_pr  libc.so.6             [.] __memset_avx512_unaligned_erms
0.18%  bench_native_pr  bench_native_prof     [.] void std::__merge_adaptive<__gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, long, int*, __gnu_cxx::__ops::_Iter_comp_iter<gkpbp::BpOsd::osd(unsigned char const*, double const*, unsigned char*, gkpbp::Result&)::{lambda(int, int)#1}> >(__gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, __gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, __gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<int> > >, long, long, int*, __gnu_cxx::__ops::_Iter_comp_iter<gkpbp::BpOsd::osd(unsigned char const*, double const*, unsigned char*, gkpbp::Result&)::{lambda(int, int)#1}>)
0.04%  bench_native_pr  libc.so.6             [.] __memmove_avx512_unaligned_erms
0.02%  bench_native_pr  bench_native_prof     [.] memmove@plt
```
`perf stat`:
```
WARNING: perf not found for kernel 6.18.33.2-microsoft

  You may need to install the following packages for this specific kernel:
    linux-tools-6.18.33.2-microsoft-standard-WSL2
    linux-cloud-tools-6.18.33.2-microsoft-standard-WSL2

  You may also want to install one of the following packages to keep up to date:
    linux-tools-standard-WSL2
    linux-cloud-tools-standard-WSL2
```

## callgrind: functions by instruction count (portable build, N = 400, 50 shots)
```
--------------------------------------------------------------------------------
Profile data file 'docs/profiling/callgrind.out' (creator: callgrind-3.22.0)
--------------------------------------------------------------------------------
I1 cache: 
D1 cache: 
LL cache: 
Timerange: Basic block 0 - 69414413
Trigger: Program termination
Profiled target:  ./build/bench_native_portable build/inst16_small.bin 1 1 (PID 1614, part 1)
Events recorded:  Ir
Events shown:     Ir
Event sort order: Ir
Thresholds:       99
Include dirs:     
User annotated:   
Auto-annotation:  on

--------------------------------------------------------------------------------
Ir                   
--------------------------------------------------------------------------------
411,289,846 (100.0%)  PROGRAM TOTALS

--------------------------------------------------------------------------------
Ir                    file:function
--------------------------------------------------------------------------------
323,244,029 (78.59%)  cpp/include/gkpbp/decoder.hpp:gkpbp::BpOsd::decode(unsigned char const*, double const*, unsigned char*)
 44,145,938 (10.73%)  /usr/include/c++/13/bits/stl_vector.h:gkpbp::BpOsd::decode(unsigned char const*, double const*, unsigned char*)
 27,928,638 ( 6.79%)  cpp/include/gkpbp/decoder.hpp:gkpbp::BpOsd::osd(unsigned char const*, double const*, unsigned char*, gkpbp::Result&)
  9,698,801 ( 2.36%)  cpp/include/gkpbp/decoder.hpp:gkpbp::BpOsd::osd_cs(double const*, unsigned char*, int)
    758,727 ( 0.18%)  ./string/../sysdeps/x86_64/multiarch/memset-vec-unaligned-erms.S:__memset_avx2_unaligned_erms
```
Hottest source lines (more than 8M instructions):
```
411,289,846 (100.0%)  PROGRAM TOTALS
323,244,029 (78.59%)  cpp/include/gkpbp/decoder.hpp:gkpbp::BpOsd::decode(unsigned char const*, double const*, unsigned char*)
44,145,938 (10.73%)  /usr/include/c++/13/bits/stl_vector.h:gkpbp::BpOsd::decode(unsigned char const*, double const*, unsigned char*)
27,928,638 ( 6.79%)  cpp/include/gkpbp/decoder.hpp:gkpbp::BpOsd::osd(unsigned char const*, double const*, unsigned char*, gkpbp::Result&)
9,698,801 ( 2.36%)  cpp/include/gkpbp/decoder.hpp:gkpbp::BpOsd::osd_cs(double const*, unsigned char*, int)
21,695,808 ( 5.28%)                  for (int e = b; e < f; ++e) {
26,427,027 ( 6.43%)                      if (v <= 0) neg ^= 1;
26,412,288 ( 6.42%)                      if (a < min1) { min2 = min1; min1 = a; arg = e; } else if (a < min2) { min2 = a; }
19,809,216 ( 4.82%)                  for (int e = b; e < f; ++e) {
19,809,216 ( 4.82%)                      const int s = neg ^ (v2c_[e] <= 0 ? 1 : 0);
21,695,808 ( 5.28%)                      const double mag = (e == arg) ? min2 : min1;
40,971,304 ( 9.96%)                      c2v_[e] = (s ? -1.0 : 1.0) * alpha_ * mag;
38,930,612 ( 9.47%)                  for (int k = col_ptr_[j]; k < col_ptr_[j + 1]; ++k) t += c2v_[col_edges_[k]];
56,047,504 (13.63%)                  for (int k = col_ptr_[j]; k < col_ptr_[j + 1]; ++k) { const int e = col_edges_[k]; v2c_[e] = t - c2v_[e]; }
41,211,005 (10.02%)  => cpp/include/gkpbp/decoder.hpp:gkpbp::BpOsd::osd(unsigned char const*, double const*, unsigned char*, gkpbp::Result&) (49x)
9,681,906 ( 2.35%)                  if (rr[w] & bit) for (int k = w; k < W_; ++k) rr[k] ^= pr[k];   // pivot row is zero left of column p
9,892,739 ( 2.41%)  => cpp/include/gkpbp/decoder.hpp:gkpbp::BpOsd::osd_cs(double const*, unsigned char*, int) (49x)
39,453,278 ( 9.59%)  	return *(this->_M_impl._M_start + __n);
406,789,267 (98.91%)  events annotated
```

## callgrind at N = 1600 (8 shots): how the split changes with size
```
Ir                    file:function
--------------------------------------------------------------------------------
207,518,862 (59.02%)  cpp/include/gkpbp/decoder.hpp:gkpbp::BpOsd::decode(unsigned char const*, double const*, unsigned char*)
 88,060,243 (25.05%)  cpp/include/gkpbp/decoder.hpp:gkpbp::BpOsd::osd(unsigned char const*, double const*, unsigned char*, gkpbp::Result&)
 28,305,304 ( 8.05%)  /usr/include/c++/13/bits/stl_vector.h:gkpbp::BpOsd::decode(unsigned char const*, double const*, unsigned char*)
 21,015,995 ( 5.98%)  cpp/include/gkpbp/decoder.hpp:gkpbp::BpOsd::osd_cs(double const*, unsigned char*, int)
    735,064 ( 0.21%)  ./string/../sysdeps/x86_64/multiarch/memset-vec-unaligned-erms.S:__memset_avx2_unaligned_erms
    668,250 ( 0.19%)  /usr/include/c++/13/bits/stl_algobase.h:gkpbp::BpOsd::osd(unsigned char const*, double const*, unsigned char*, gkpbp::Result&)
    631,800 ( 0.18%)  /usr/include/c++/13/bits/move.h:gkpbp::BpOsd::osd(unsigned char const*, double const*, unsigned char*, gkpbp::Result&)
    490,784 ( 0.14%)  ./elf/./elf/dl-lookup.c:do_lookup_x
    485,870 ( 0.14%)  ./elf/../sysdeps/generic/dl-new-hash.h:_dl_lookup_symbol_x
    432,756 ( 0.12%)  /usr/include/c++/13/bits/stl_algo.h:int* std::__move_merge<__gnu_cxx::__normal_iterator<int*, std::vector<int, std::allocator<in
```

## memcheck
```
==1653==     in use at exit: 0 bytes in 0 blocks
==1653==   total heap usage: 327 allocs, 327 frees, 530,036 bytes allocated
==1653== All heap blocks were freed -- no leaks are possible
==1653== ERROR SUMMARY: 0 errors from 0 contexts (suppressed: 0 from 0)
```

## cachegrind
```
I refs:        414,536,213
I1  misses:          2,673
LLi misses:          2,494
D refs:        104,928,085  (85,795,930 rd   + 19,132,155 wr)
D1  misses:      1,320,343  (   934,054 rd   +    386,289 wr)
LLd misses:         13,081  (     7,586 rd   +      5,495 wr)
D1  miss rate:         1.3% (       1.1%     +        2.0%  )
LL misses:          15,575  (    10,080 rd   +      5,495 wr)
LL miss rate:          0.0% (       0.0%     +        0.0%  )
```

## After the change
Alternating old/new timings (N = 1600, one thread, five pairs, microseconds per shot):
```
old: 4138.5 us/shot
new: 3261.0 us/shot
old: 4164.9 us/shot
new: 3249.0 us/shot
old: 4123.8 us/shot
new: 3234.3 us/shot
old: 4127.0 us/shot
new: 3237.8 us/shot
old: 4144.7 us/shot
new: 3274.6 us/shot
```

## Reading the profile
Two complementary views. Callgrind instruction counts are deterministic for a given compiler (gcc 13.3, portable `-O2` build, hard prior, sigma = 0.50). `perf` is wall-clock sampling on the optimised `-O3 -march=native` build (N = 1600, 300 shots, one thread).

- **Belief propagation dominates at small N.** At N = 400 about 89 % of instructions are in the BP loops (78.6 % in `decode` plus 10.7 % attributed to inlined `std::vector` indexing) and about 9 % in OSD. BP converges on only about 5 % of these shots, so nearly every shot runs all 100 iterations before OSD is reached.
- **OSD grows with code size.** At N = 1600 callgrind gives about 67 % BP and about 31 % OSD by instruction count. By time, `perf` puts OSD at about half: BP (inlined into the OpenMP region) 46.8 %, `osd_cs` 28.6 %, `osd` 22.4 %. So OSD costs more per instruction than BP. My untested guess is bit-matrix elimination and memory access on one side, and the vectorised optimised build on the other; I have not checked this.
- **Inside BP the cost is spread over the indirect-index loops** of the two node updates: the variable-node write-back `v2c_[e] = t - c2v_[e]` (13.6 %), the check-node message write `c2v_[e] = ...` (10.0 %), the posterior sum `t += c2v_[col_edges_[k]]` (9.5 %), and the min search (about 13 %). All of them walk edges through `row_ptr_` or `col_edges_`.
- **It is not a last-level-cache problem at this size.** Cachegrind simulates the host's caches. On the author's laptop (WSL2) it reports a 1.3 % D1 miss rate and about 0 % last-level misses at N = 400 (on a machine with a larger L1 the D1 rate was 0.1 %), so the working set stays in cache and the time is mostly instruction count and indirection rather than memory stalls. I have not measured cache behaviour at N = 1600.
- **Memory safety:** memcheck reports 0 errors and no leaks (327 allocations, all freed).

## Change made
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

## Caveats
Valgrind runs a portable `-O2` build without `-march=native` and without OpenMP (valgrind cannot decode some AVX-512 instructions and OpenMP produces noise); timings use the optimised build.

#!/usr/bin/env bash
# One-command profiling run for the C++ decoder core (Linux / WSL). Run from the repo root inside the venv:
#     bash scripts/profile_all.sh
# Optional: PERF=/usr/lib/linux-tools/<version>/perf bash scripts/profile_all.sh   (WSL often needs the full path)
# Writes everything to docs/profiling/. Steps that cannot run on this machine (e.g. perf under WSL) are skipped with a note.
set -u
OUT=docs/profiling; mkdir -p "$OUT" build
PERF=${PERF:-perf}
INC="-Icpp/include"; SRC="cpp/bench/bench_native.cpp"
note() { echo "== $*"; }

note "machine"; { uname -a; lscpu | grep -E "Model name|^CPU\(s\)|Thread|L1d|L2|L3"; g++ --version | head -1; python --version; } | tee "$OUT/machine.txt"

note "build: optimised+symbols (for timing/perf) and portable (for valgrind)"
g++ -O3 -g -fno-omit-frame-pointer -std=c++17 -fopenmp -march=native -DNDEBUG $INC $SRC -o build/bench_native_prof || exit 1
# valgrind cannot decode every -march=native instruction (AVX-512 -> "Illegal instruction") and reports OpenMP runtime noise,
# so the valgrind binary is built without both:
g++ -O2 -g -std=c++17 $INC $SRC -o build/bench_native_portable || exit 1

note "inputs"
python scripts/export_instance.py 24 0.50 300 hard 7 build/inst24.bin
python scripts/export_instance.py 32 0.50 300 hard 7 build/inst32.bin
python scripts/export_instance.py 16 0.50 50  hard 7 build/inst16_small.bin

note "timing (1 thread, best of 3)"
{ ./build/bench_native_prof build/inst24.bin 1 3; ./build/bench_native_prof build/inst32.bin 1 3; } | tee "$OUT/timing.txt"

note "perf (hardware counters often unavailable under WSL; falls back to software clock)"
if command -v "$PERF" >/dev/null 2>&1 || [ -x "$PERF" ]; then
  "$PERF" stat -d ./build/bench_native_prof build/inst24.bin 1 3 > "$OUT/perf_stat.txt" 2>&1 || echo "perf stat failed (see file)"
  if "$PERF" record -g -o "$OUT/perf.data" ./build/bench_native_prof build/inst32.bin 1 3 >/dev/null 2>&1 ||
     "$PERF" record -e cpu-clock -g -o "$OUT/perf.data" ./build/bench_native_prof build/inst32.bin 1 3 >/dev/null 2>&1; then
    "$PERF" report -i "$OUT/perf.data" --stdio --no-children 2>/dev/null | head -80 > "$OUT/perf_report.txt"; echo "perf report saved"
  else echo "perf record failed -> use the callgrind results as the profiler"; fi
else echo "perf not installed -> skipping (callgrind below is the profiler)"; fi

note "callgrind (instruction-level hotspots)"
if command -v valgrind >/dev/null 2>&1; then
  valgrind --tool=callgrind --callgrind-out-file="$OUT/callgrind.out" ./build/bench_native_portable build/inst16_small.bin 1 1 > "$OUT/callgrind_run.txt" 2>&1
  callgrind_annotate "$OUT/callgrind.out" 2>/dev/null | sed -E 's/ \[[^]]*\]//g' > "$OUT/callgrind_functions.txt"
  callgrind_annotate --auto=yes "$OUT/callgrind.out" 2>/dev/null | awk '/^ *[0-9,]+ \(/ {t=$1; gsub(",","",t); if (t+0 > 8000000) print}' | sed 's/^ *//' | sed -E 's/ \[[^]]*\]//g' > "$OUT/callgrind_hot_lines.txt"
  note "callgrind at N = 1600 (shows how the BP / OSD split changes with size)"
  python scripts/export_instance.py 32 0.50 8 hard 7 build/inst32_small.bin >/dev/null
  valgrind --tool=callgrind --callgrind-out-file="$OUT/callgrind_n1600.out" ./build/bench_native_portable build/inst32_small.bin 1 1 > "$OUT/callgrind_n1600_run.txt" 2>&1
  callgrind_annotate "$OUT/callgrind_n1600.out" 2>/dev/null | sed -n '/file:function/,/^$/p' | sed -E 's/ \[[^]]*\]//g' | cut -c1-150 | head -12 > "$OUT/callgrind_n1600_functions.txt"
  note "memcheck (memory errors / leaks)"
  valgrind --leak-check=full ./build/bench_native_portable build/inst16_small.bin 1 1 > "$OUT/memcheck.txt" 2>&1; grep -E "ERROR SUMMARY|in use at exit|total heap usage|All heap blocks" "$OUT/memcheck.txt"
  note "cachegrind (cache misses)"
  valgrind --tool=cachegrind --cache-sim=yes --cachegrind-out-file="$OUT/cachegrind.out" ./build/bench_native_portable build/inst16_small.bin 1 1 > "$OUT/cachegrind_run.txt" 2>&1
  grep -E "I +refs|I1  misses|LLi misses|D +refs|D1  misses|LLd misses|LL misses|LL miss rate|D1  miss rate" "$OUT/cachegrind_run.txt" | sed 's/==[0-9]*== //' | tee "$OUT/cachegrind_summary.txt"
else echo "valgrind not installed: sudo apt install valgrind"; fi
echo; echo "done. Next: python scripts/make_profiling_md.py   (builds docs/PROFILING.md from these files)"

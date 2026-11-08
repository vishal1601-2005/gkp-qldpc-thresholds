"""Assemble docs/PROFILING.md from the files written by scripts/profile_all.sh (docs/profiling/*).
Write your prose in docs/profiling/reading.md and docs/profiling/change.md; they are inserted automatically on every run."""
import pathlib
d = pathlib.Path("docs/profiling")
def rd(name, n=None):
    p = d / name
    if not p.exists(): return f"(not available on this machine: {name})"
    t = p.read_text(errors="replace").rstrip().splitlines()
    return "\n".join(t[:n] if n else t)
def prose(name, placeholder):
    p = d / name
    return p.read_text().strip() if p.exists() else placeholder
import re
_pr = rd('perf_report.txt')
_flat = [l.strip() for l in _pr.splitlines() if re.match(r'\s+[0-9.]+%', l)][:10]
perf_flat = "\n".join(_flat) if _flat else _pr
after = ""
if (d / "timing_before_after.txt").exists():
    after = f"""## After the change
Alternating old/new timings (N = 1600, one thread, five pairs, microseconds per shot):
```
{rd('timing_before_after.txt')}
```
"""
md = f"""# Profiling the C++ decoder core

Generated from `scripts/profile_all.sh`; raw files are in `docs/profiling/`.

## Machine
```
{rd('machine.txt')}
```

## Timing (single thread, best of 3; N = 900 and N = 1600, hard sigma = 0.50, OSD-CS(7), 300 shots)
```
{rd('timing.txt')}
```

## perf
```
{perf_flat}
```
`perf stat`:
```
{rd('perf_stat.txt', 40)}
```

## callgrind: functions by instruction count (portable build, N = 400, 50 shots)
```
{rd('callgrind_functions.txt', 30)}
```
Hottest source lines (more than 8M instructions):
```
{rd('callgrind_hot_lines.txt', 25)}
```

## callgrind at N = 1600 (8 shots): how the split changes with size
```
{rd('callgrind_n1600_functions.txt', 12)}
```

## memcheck
```
{chr(10).join(l for l in rd('memcheck.txt').splitlines() if any(k in l for k in ('ERROR SUMMARY','in use at exit','total heap usage','All heap blocks','definitely lost','Invalid')))}
```

## cachegrind
```
{rd('cachegrind_summary.txt')}
```

{after}
## Reading the profile
{prose('reading.md', '*(Write 3-5 sentences in docs/profiling/reading.md: which stage dominates, how it changes with N, what the cache numbers say.)*')}

## Change made
{prose('change.md', '*(Write docs/profiling/change.md: hotspot, change, before/after numbers, and the pytest result. Skip if no safe improvement was found.)*')}

## Caveats
Valgrind runs a portable `-O2` build without `-march=native` and without OpenMP (valgrind cannot decode some AVX-512 instructions and OpenMP produces noise); timings use the optimised build.
"""
pathlib.Path("docs/PROFILING.md").write_text(md); print("wrote docs/PROFILING.md")

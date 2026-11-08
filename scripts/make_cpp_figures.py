"""Figures for the C++ core and the BB cross-check: docs/fig4_decoder_speedup.png, docs/fig5_bb_code_capacity.png."""
import csv, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt

rows = list(csv.DictReader(open("data/bench/bench_ubuntu_12core_summary.csv")))
fig, ax = plt.subplots(figsize=(6.4, 3.8)); labels = [f"N={r['N']} {r['mode']} σ={r['sigma']}" for r in rows]
x = np.arange(len(rows)); w = 0.38
ld = [float(r["ldpc_ms_per_shot"]) for r in rows]; cp = [float(r["cpp_1thread_ms"]) for r in rows]
ax.bar(x - w / 2, ld, w, label="ldpc BpOsdDecoder (Python loop)"); ax.bar(x + w / 2, cp, w, label="C++ core, 1 thread")
ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=6, rotation=55, ha='right'); ax.set_ylabel("ms per shot (log)")
ax.set_title("BP + OSD-CS(7), identical shots/priors, single core", fontsize=9); ax.legend(fontsize=7)
for xi, a_, b_ in zip(x, ld, cp): ax.text(xi + w / 2, b_ * 1.15, f"{a_/b_:.1f}×", ha="center", fontsize=6)
fig.tight_layout(); fig.savefig("docs/fig4_decoder_speedup.png", dpi=170)

A = np.loadtxt("data/raw/bb_code_capacity_results.csv", delimiter=",")
fig, axs = plt.subplots(1, 2, figsize=(9, 3.6), sharey=True)
for ax, mode, name in zip(axs, (0, 1), ("hard", "soft")):
    for N, c in ((72, "C0"), (90, "C1"), (108, "C2"), (144, "C3"), (288, "C4")):
        M = A[(A[:, 0] == mode) & (A[:, 1] == N)]; M = M[np.argsort(M[:, 3])]
        k, d = int(M[0, 7]), int(M[0, 8]); ax.plot(M[:, 3], M[:, 5] / M[:, 4], "o-", ms=3, color=c, label=f"[[{N},{k},{d}]]")
    ax.set_title(f"BB codes, {name}-decision BP+OSD-CS(7)"); ax.set_xlabel("GKP shift std σ"); ax.grid(alpha=.3)
axs[0].set_ylabel("block failure probability"); axs[0].legend(fontsize=7); fig.tight_layout(); fig.savefig("docs/fig5_bb_code_capacity.png", dpi=170)

# fig6: multi-core results measured on a 12-logical-core machine (left: latest bench_decoders run; right: native scaling file)
import re
rows12 = list(csv.DictReader(open("data/bench/bench_ubuntu_12core_summary.csv")))
fig, axs = plt.subplots(1, 2, figsize=(9, 3.6))
for N, c in (("400", "C0"), ("900", "C1"), ("1600", "C2")):
    r = [x for x in rows12 if x["N"] == N and x["mode"] == "hard" and x["sigma"] == "0.50"][0]
    ld = float(r["ldpc_ms_per_shot"]); axs[0].plot([1, 4, 8], [ld / float(r[k]) for k in ("cpp_1thread_ms", "cpp_4threads_ms", "cpp_8threads_ms")], "o-", color=c, label=f"N={N}")
axs[0].set_xlabel("C++ threads"); axs[0].set_ylabel("speedup over ldpc Python loop"); axs[0].set_title("hard σ=0.50, BP+OSD-CS(7)", fontsize=9); axs[0].legend(); axs[0].grid(alpha=.3)
us = {}
for line in open("data/bench/native_cpp_scaling_12core.txt"):
    m = re.search(r"threads=(\d+).*?([0-9.]+) us/shot", line)
    if m: us[int(m.group(1))] = float(m.group(2))
ts = sorted(us); axs[1].plot(ts, [us[1] / us[t] for t in ts], "o-", label="measured (N=1600, native C++)"); axs[1].plot(ts, ts, "k--", lw=.8, label="ideal")
axs[1].set_xlabel("threads"); axs[1].set_ylabel("speedup vs 1 thread"); axs[1].set_title("native C++ scaling (earlier run)", fontsize=9); axs[1].legend(fontsize=7); axs[1].grid(alpha=.3)
fig.tight_layout(); fig.savefig("docs/fig6_multicore_scaling.png", dpi=170)
print("figures written")

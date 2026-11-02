"""Compact view of layer_a_results.csv: block failure P (aggregated over codes) with the min-max across codes, plus per-code distances.
Run: python show_table_v2.py      then paste the whole output into the chat."""
import numpy as np
A = np.loadtxt("layer_a_results.csv", delimiter=",", ndmin=2)
for mode, name in ((0, "HARD"), (1, "SOFT")):
    M = A[A[:, 0] == mode]
    if not len(M): continue
    ns = sorted(set(M[:, 1].astype(int))); sg = sorted(set(np.round(M[:, 3], 4)))
    print(f"\n===== {name}: block failure P  [min-max over codes] =====")
    print("codes (n:seed=d1): " + "  ".join(f"{int(n)}:{int(s)}=d{int(d)}" for n, s, d in sorted({(r[1], r[2], r[8]) for r in M})))
    print("sigma   " + "".join(f"{'n='+str(n):>20}" for n in ns))
    for s in sg:
        row = f"{s:.3f}  "
        for n in ns:
            r = M[(M[:, 1] == n) & np.isclose(M[:, 3], s)]
            if not len(r): row += f"{'-':>20}"; continue
            P = r[:, 5] / r[:, 4]
            row += f"{r[:,5].sum()/r[:,4].sum():>9.3f} [{P.min():.2f}-{P.max():.2f}]"
        print(row)

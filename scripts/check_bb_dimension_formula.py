"""Usage: python scripts/check_bb_dimension_formula.py     (needs: pip install galois; about a minute)

Checks results of Sabo, Can, Marquis, "Spectral Theory of Semisimple Bivariate Bicycle Codes" (arXiv:2608.27565)
against this repo's bivariate_bicycle() and plain GF(2) linear algebra:

 A. Theorem 81: for a semisimple ring (l, m odd, q = 2), k = 2 * |{(i, j) in Z_l x Z_m : a(alpha^i, beta^j) = b(alpha^i, beta^j) = 0}|,
    where alpha, beta are primitive l-th and m-th roots of unity in an extension of F_2. Random codes, coprime or not.
 B. The explicit codes of the paper's Examples 96, 132, 133 and 134 (stated k = 2, 8, 24, 10).
 C. Prop. 65 (coprime l, m, any parity): k = 2 * deg gcd(a(z), b(z), z^(lm) - 1) with x -> z^m, y -> z^l.
"""
import random
import numpy as np
import galois
from gkp_qldpc.codes import bivariate_bicycle

def gf2_rank(M):
    rows = [int("".join(map(str, r)), 2) for r in (np.asarray(M) % 2).astype(int)]
    rank = 0
    while rows:
        p = rows.pop()
        if p == 0: continue
        rank += 1
        low = p & -p
        rows = [r ^ p if r & low else r for r in rows]
    return rank

def k_numeric(l, m, a, b):
    HX, HZ = bivariate_bicycle(l, m, a, b)
    return HX.shape[1] - gf2_rank(HX.toarray()) - gf2_rank(HZ.toarray())

_fields = {}
def field_for(l, m):
    """Smallest GF(2^e) with l | 2^e - 1 and m | 2^e - 1 (l, m odd)."""
    e = 1
    while (2 ** e - 1) % l or (2 ** e - 1) % m: e += 1
    if e not in _fields: _fields[e] = galois.GF(2 ** e)
    return _fields[e], e

def common_zeros(l, m, a, b):
    GF, e = field_for(l, m)
    N = 2 ** e - 1
    g = GF.primitive_element
    I, J = np.meshgrid(np.arange(l), np.arange(m), indexing="ij")
    def evaluate(terms):
        total = GF.Zeros(I.shape)
        for s, t in terms:               # alpha^(i s) * beta^(j t) = g^((N/l)(i s) + (N/m)(j t))
            k = ((N // l) * ((I * s) % l) + (N // m) * ((J * t) % m)) % N
            total = total + g ** k
        return total
    return int(np.count_nonzero((evaluate(a) == 0) & (evaluate(b) == 0)))

# ---- Prop. 65 helpers (polynomials over F2 as int bitmasks)
def pmod(p, q):
    dq = q.bit_length()
    while p.bit_length() >= dq: p ^= q << (p.bit_length() - dq)
    return p
def pgcd(p, q):
    while q: p, q = q, pmod(p, q)
    return p
def to_z(terms, l, m):
    p = 0
    for i, j in terms: p ^= 1 << ((m * i + l * j) % (l * m))
    return p
def prop65_k(l, m, a, b):
    g = pgcd(pgcd(to_z(a, l, m), to_z(b, l, m)), (1 << (l * m)) | 1)
    return 2 * (g.bit_length() - 1) if g else None

def random_code(l, m, rng):
    a = [(rng.randrange(l), rng.randrange(m)) for _ in range(rng.choice([2, 3, 3, 4]))]
    b = [(rng.randrange(l), rng.randrange(m)) for _ in range(rng.choice([2, 3, 3, 4]))]
    return a, b

rng = random.Random(1)

print("A. Theorem 81 on random semisimple codes (l, m odd; k = 2 x number of common zeros)")
pairs = [(3, 3), (3, 5), (3, 7), (5, 5), (3, 9), (3, 15), (5, 9), (7, 7), (9, 9), (5, 15), (7, 9), (3, 11)]
ok = tot = 0; ks = []
for l, m in pairs:
    for _ in range(12):
        a, b = random_code(l, m, rng)
        kn, kp = k_numeric(l, m, a, b), 2 * common_zeros(l, m, a, b)
        tot += 1; ok += (kn == kp); ks.append(kn)
        if kn != kp: print("   MISMATCH", l, m, a, b, kn, kp)
print(f"   {ok}/{tot} match; k ranged over {sorted(set(ks))}")

print("B. Explicit codes from the paper")
examples = {
    "Example 96  [[42, 2]]":  (3, 7,  [(1, 1), (1, 4), (2, 1), (2, 3)], [(0, 2), (2, 3), (2, 4), (2, 6)], 2),
    "Example 132 [[90, 8]]":  (3, 15, [(0, 0), (0, 1), (0, 5)], [(2, 0), (1, 0), (0, 3)], 8),
    "Example 133 [[162, 24]]": (9, 9, [(0, 0), (0, 1), (0, 2)], [(0, 3), (3, 0), (6, 0)], 24),
    "Example 134 [[434, 10]]": (31, 7, [(0, 0), (6, 1), (27, 4)], [(0, 0), (15, 6), (24, 3)], 10),
}
for name, (l, m, a, b, k_paper) in examples.items():
    kz = 2 * common_zeros(l, m, a, b)
    kn = k_numeric(l, m, a, b) if l * m <= 300 else None   # the 434-qubit rank is fine too, but keep runtime small
    print(f"   {name}: paper k = {k_paper}; Theorem 81 count = {kz}; GF(2) rank = {kn if kn is not None else 'skipped'}")

print("C. Prop. 65 on random coprime (l, m), any parity")
ok = tot = 0
for l, m in [(3, 5), (3, 7), (5, 7), (3, 11), (5, 9), (7, 9), (4, 5), (2, 7), (4, 9), (6, 5)]:
    for _ in range(20):
        a, b = random_code(l, m, rng)
        if to_z(a, l, m) == 0 and to_z(b, l, m) == 0: continue
        tot += 1; ok += (k_numeric(l, m, a, b) == prop65_k(l, m, a, b))
print(f"   {ok}/{tot} match")

"""Code constructions (hypergraph product, toric), GF(2) helpers, and the screened random (3,4)-LDPC ensemble.

The screened ensemble exists because finite-size threshold estimates for constant-rate qLDPC families are dominated by a
"code lottery": random (3,4) codes at n<=48 range from classical distance 2 to 12, and a single d=2 code ruins a crossing.
make_good_code() therefore requires full row rank (fixed k) and classical distance >= 2*ceil(n/8) (see dmin_rule).
"""
import itertools
import numpy as np
import scipy.sparse as sp
import ldpc.mod2 as m2

def rand_regular(n, dv, dc, rng, max_tries=100000):
    """Random (dv,dc)-regular parity-check matrix H (m x n), m = n*dv/dc, no repeated edges."""
    assert (n * dv) % dc == 0
    m = n * dv // dc
    for _ in range(max_tries):
        stubs = np.repeat(np.arange(n), dv)
        rng.shuffle(stubs)
        rows = np.arange(len(stubs)) // dc
        pairs = rows * n + stubs
        if len(np.unique(pairs)) != len(pairs):
            continue
        H = np.zeros((m, n), dtype=np.uint8)
        H[rows, stubs] = 1
        return H
    raise RuntimeError("no simple regular graph found")

def hgp(H1, H2):
    """Hypergraph product (Tillich-Zemor). Returns (HX, HZ) as scipy CSR uint8, N = n1 n2 + m1 m2."""
    H1 = sp.csr_matrix(H1, dtype=np.int64); H2 = sp.csr_matrix(H2, dtype=np.int64)
    m1, n1 = H1.shape; m2_, n2 = H2.shape
    HX = sp.hstack([sp.kron(H1, sp.identity(n2)), sp.kron(sp.identity(m1), H2.T)]).tocsr()
    HZ = sp.hstack([sp.kron(sp.identity(n1), H2), sp.kron(H1.T, sp.identity(m2_))]).tocsr()
    HX.data %= 2; HZ.data %= 2
    HX.eliminate_zeros(); HZ.eliminate_zeros()
    assert ((HX @ HZ.T).toarray() % 2 == 0).all(), "HX HZ^T != 0"
    return HX.astype(np.uint8), HZ.astype(np.uint8)

def ring_rep(L):
    """Parity checks of the periodic repetition code (L x L)."""
    H = np.zeros((L, L), dtype=np.uint8)
    for i in range(L):
        H[i, i] = 1; H[i, (i + 1) % L] = 1
    return H

def toric(L):
    return hgp(ring_rep(L), ring_rep(L))

class LogicalChecker:
    """For X-type errors detected by HZ: residual r (HZ r = 0) is a stabilizer iff K r = 0,
    where K spans ker(HX) (r in rowspace(HX)  <=>  r orthogonal to ker(HX))."""
    def __init__(self, HX):
        K = m2.nullspace(sp.csr_matrix(HX))
        self.K = sp.csr_matrix(K).astype(np.uint8)
    def is_logical(self, r):
        return bool(((self.K @ r) % 2).any())

def code_params(HX, HZ):
    N = HX.shape[1]
    k = N - m2.rank(sp.csr_matrix(HX)) - m2.rank(sp.csr_matrix(HZ))
    return N, int(k)

def kernel_dense(HX):
    """Dense float32 basis of ker(HX) over GF(2) (used for batched logical-failure checks)."""
    K = m2.nullspace(sp.csr_matrix(HX))
    return np.asarray(K.todense(), dtype=np.float32)

def logical_fail_batch(R, Kd):
    """R: (B,N) residual X-errors with HZ r = 0.  Returns bool array: True if r is a nontrivial logical."""
    return (((R.astype(np.float32) @ Kd.T) % 2) > 0.5).any(axis=1)

def make_code(n, seed, dv=3, dc=4):
    rng = np.random.default_rng(seed)
    H = rand_regular(n, dv, dc, rng)
    HX, HZ = hgp(H, H)
    return HX, HZ

def dmin_rule(n):
    """Minimum classical distance required of a sampled code: d1 >= 2*ceil(n/8) (4,4,6,6,8,8,10,10,12,12 for n=12,16,...,48).
    Identical to the earlier step rule for n<=28, so existing results stay valid. Rejects finite-size 'bad lottery' codes."""
    return 2 * int(np.ceil(n / 8))

def classical_distance(H):
    K = np.asarray(m2.nullspace(sp.csr_matrix(H)).todense(), dtype=np.uint8)
    C = np.array(list(itertools.product([0, 1], repeat=K.shape[0])), dtype=np.uint8)[1:]
    return int(((C @ K) % 2).sum(1).min())

def make_good_code(n, seed, dv=3, dc=4):
    """Random (dv,dc)-regular H with FULL row rank (so k = (n/4)^2 exactly, no transpose logicals) and d1 >= dmin_rule(n).
    Deterministic in (n, seed). Returns HX, HZ, k, d1."""
    rng = np.random.default_rng(seed)
    for _ in range(100000):
        H = rand_regular(n, dv, dc, rng)
        if m2.rank(sp.csr_matrix(H)) < H.shape[0]: continue
        d1 = classical_distance(H)
        if d1 < dmin_rule(n): continue
        HX, HZ = hgp(H, H)
        return HX, HZ, (n // 4) ** 2, d1
    raise RuntimeError("no code found")


# ---------------- bivariate-bicycle (BB) codes (Bravyi et al., Nature 627, 778 (2024)) ----------------
def bivariate_bicycle(l, m, a_terms, b_terms):
    """BB code on N = 2*l*m qubits. Polynomials A = sum x^i y^j, B = sum x^i y^j in the ring F2[x,y]/(x^l - 1, y^m - 1),
    with x = S_l (x) I_m and y = I_l (x) S_m (cyclic shifts). H_X = [A | B], H_Z = [B^T | A^T].
    a_terms / b_terms: lists of (i, j) exponents. Returns (HX, HZ) as CSR uint8."""
    def shift(k):
        return sp.csr_matrix(np.roll(np.eye(k, dtype=np.int64), 1, axis=1))
    I = sp.identity(l * m, format="csr", dtype=np.int64)
    x = sp.kron(shift(l), sp.identity(m, dtype=np.int64), format="csr")
    y = sp.kron(sp.identity(l, dtype=np.int64), shift(m), format="csr")
    def power(M, k):
        out = I
        for _ in range(k): out = out @ M
        return out
    def poly(terms):
        out = sp.csr_matrix((l * m, l * m), dtype=np.int64)
        for i, j in terms: out = out + power(x, i) @ power(y, j)
        out.data %= 2; out.eliminate_zeros(); return out
    A, B = poly(a_terms), poly(b_terms)
    HX = sp.hstack([A, B]).tocsr(); HZ = sp.hstack([B.T, A.T]).tocsr()
    assert not ((HX @ HZ.T).toarray() % 2).any(), "HX HZ^T != 0"
    return HX.astype(np.uint8), HZ.astype(np.uint8)


# [[N, k, d]] as published; N and k are verified in the tests, d is quoted from the literature (not recomputed here).
BB_CODES = {
    "72_12_6":    dict(l=6,  m=6,  a=[(3, 0), (0, 1), (0, 2)], b=[(0, 3), (1, 0), (2, 0)], d=6),
    "90_8_10":    dict(l=15, m=3,  a=[(9, 0), (0, 1), (0, 2)], b=[(0, 0), (2, 0), (7, 0)], d=10),
    "108_8_10":   dict(l=9,  m=6,  a=[(3, 0), (0, 1), (0, 2)], b=[(0, 3), (1, 0), (2, 0)], d=10),
    "144_12_12":  dict(l=12, m=6,  a=[(3, 0), (0, 1), (0, 2)], b=[(0, 3), (1, 0), (2, 0)], d=12),
    "288_12_18":  dict(l=12, m=12, a=[(3, 0), (0, 2), (0, 7)], b=[(0, 3), (1, 0), (2, 0)], d=18),
}


def bb_code(name):
    c = BB_CODES[name]; HX, HZ = bivariate_bicycle(c["l"], c["m"], c["a"], c["b"]); return HX, HZ, c["d"]

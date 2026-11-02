"""Finite-size analysis: pairwise crossings with cluster-bootstrap errors, and a collapse fit (used for the toric-code gate)."""
import numpy as np
from scipy.optimize import least_squares

def fss_fit(Ls, ps, fails, shots, nu0=1.5, seed=0, n_boot=300):
    """Fit p_L = a0 + a1 x + a2 x^2, x=(p-pc) L^{1/nu}. Returns (pc, nu, pc_err, nu_err) with parametric bootstrap.
    Ls,ps,fails,shots are equal-length 1d arrays (one entry per (L,p) point)."""
    from scipy.optimize import least_squares
    Ls = np.asarray(Ls, float); ps = np.asarray(ps, float)
    fails = np.asarray(fails, float); shots = np.asarray(shots, float)
    def fit(f):
        y = f / shots
        yv = (f + 0.5) / (shots + 1.0)                      # shrunk estimate: keeps error bars finite at 0 / all failures
        err = np.sqrt(yv * (1 - yv) / shots)
        def res(th):
            pc, nu, a0, a1, a2 = th
            x = (ps - pc) * Ls ** (1 / nu)
            return (a0 + a1 * x + a2 * x ** 2 - y) / err
        th0 = [ps.mean(), nu0, 0.3, 1.0, 0.0]
        lo = [ps.min(), 0.8, -1, -50, -50]; hi = [ps.max(), 4.0, 2, 50, 50]
        r = least_squares(res, th0, bounds=(lo, hi))
        return r.x, (r.fun ** 2).sum() / max(len(y) - 5, 1)
    th, chi2 = fit(fails)
    rng = np.random.default_rng(seed); boots = []
    for _ in range(n_boot):
        fb = rng.binomial(shots.astype(int), np.clip(fails / shots, 0, 1))
        try: boots.append(fit(fb)[0][:2])
        except Exception: pass
    boots = np.array(boots)
    if min(abs(th[0] - ps.min()), abs(th[0] - ps.max())) < 1e-6 or th[1] in (0.8, 4.0):
        print("   NOTE: fit hit a bound -> no reliable crossing inside the scanned window; widen/shift the window")
    return th[0], th[1], boots[:, 0].std(), boots[:, 1].std(), chi2

def block_curve(M, n, rng=None):
    """Aggregated block-failure curve {sigma: P} for code size n from a results array with columns
    (mode, n, codeseed, sigma, shots, fails, syndrome_mismatch, k, d1). Optionally cluster-bootstrap over codes."""
    Mn = M[M[:, 1] == n]; seeds = sorted(set(Mn[:, 2].astype(int)))
    if rng is not None: seeds = list(rng.choice(seeds, size=len(seeds), replace=True))
    out = {}
    for s in sorted(set(np.round(Mn[:, 3], 4))):
        f = sh = 0
        for sd in seeds:
            r = Mn[(Mn[:, 2] == sd) & np.isclose(Mn[:, 3], s)]
            if len(r): f += r[0, 5]; sh += r[0, 4]
        if sh:
            if rng is not None: f = rng.binomial(int(sh), min(max(f / sh, 0), 1))
            out[s] = f / sh
    return out


def crossing(ca, cb):
    """sigma where the larger-code curve cb rises above the smaller-code curve ca (linear interpolation); nan if none."""
    common = sorted(set(ca) & set(cb))
    d = np.array([cb[s] - ca[s] for s in common]); sg = np.array(common)
    for i in range(len(d) - 1):
        if d[i] < 0 <= d[i + 1]: return sg[i] + (sg[i + 1] - sg[i]) * (-d[i]) / (d[i + 1] - d[i])
    return np.nan


def classify_pair(dist_a, dist_b):
    """'clean' (one distance per size, increasing), 'same' (equal distance: multiplicity-dominated), or 'mixed'."""
    if len(dist_a) > 1 or len(dist_b) > 1: return "mixed"
    return "same" if dist_a[0] == dist_b[0] else "clean"


def pairwise_crossings(M, n_boot=300, seed=7):
    """For a results array M (one decoder mode) return a list of dicts: n_a, n_b, sigma_x, err, n_found, category."""
    ns = sorted(set(M[:, 1].astype(int)))
    dist = {n: sorted({int(x) for x in M[M[:, 1] == n][:, 8]}) for n in ns}
    res = []
    for a, b in zip(ns[:-1], ns[1:]):
        c0 = crossing(block_curve(M, a), block_curve(M, b))
        rng = np.random.default_rng(seed)
        bs = [crossing(block_curve(M, a, rng), block_curve(M, b, rng)) for _ in range(n_boot)]
        ok = [x for x in bs if np.isfinite(x)]
        res.append(dict(n_a=a, n_b=b, sigma_x=float(c0), err=float(np.std(ok)) if len(ok) > 30 else float("nan"),
                        n_found=len(ok), category=classify_pair(dist[a], dist[b]), dist_a=dist[a], dist_b=dist[b]))
    return res

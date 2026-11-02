"""First-principles check of the GKP loss channel (finite-energy GKP state, exact Kraus loss, exact homodyne marginal).

The x-marginal after pure loss is exactly the convolution of the scaled input marginal with the vacuum marginal, so this
verifies the effective-variance model of the *bit-flip probability*. It does NOT characterise the full logical channel:
pure loss on approximate GKP qubits is not a Pauli channel (Harris et al., PRA 112, 042417).
"""
import numpy as np
from scipy.special import gammaln
from .channel import SQRTPI, p_flip_exact, sigma2_eff, sigma_from_p

def hermite_funcs(x, N):
    """psi_n(x), n=0..N-1, stable forward recurrence. Returns (len(x), N)."""
    psi = np.zeros((len(x), N))
    psi[:, 0] = np.pi ** -0.25 * np.exp(-x ** 2 / 2)
    if N > 1: psi[:, 1] = np.sqrt(2.0) * x * psi[:, 0]
    for n in range(1, N - 1):
        psi[:, n + 1] = np.sqrt(2.0 / (n + 1)) * x * psi[:, n] - np.sqrt(n / (n + 1.0)) * psi[:, n - 1]
    return psi

def gkp_fock(Delta2, N, S=60):
    """|0_D> = e^{-D^2 n} sum_s |x = 2 s sqrt(pi)>  in Fock basis (real coefficients), normalized."""
    xs = 2 * SQRTPI * np.arange(-S, S + 1)
    c = hermite_funcs(xs, N).sum(0) * np.exp(-Delta2 * np.arange(N))
    return c / np.linalg.norm(c)

def apply_loss(rho, eta):
    """Exact pure-loss channel on a (truncated) density matrix, Kraus form."""
    N = rho.shape[0]
    out = np.zeros_like(rho)
    m = np.arange(N)
    for k in range(N):
        dim = N - k
        # A_k |n> = sqrt(C(n,k)) eta^{(n-k)/2} (1-eta)^{k/2} |n-k>
        n = m[k:]
        a = np.exp(0.5 * (gammaln(n + 1) - gammaln(k + 1) - gammaln(n - k + 1))
                   + 0.5 * (n - k) * np.log(eta) + 0.5 * k * np.log1p(-eta))
        out[:dim, :dim] += (a[:, None] * a[None, :]) * rho[k:, k:]
    return out

def xmarginal(rho, X, dx):
    x = np.arange(-X, X + dx / 2, dx)
    psi = hermite_funcs(x, rho.shape[0])
    P = np.einsum('ij,jk,ik->i', psi, rho, psi, optimize=True)
    return x, P


def logical_error_vs_loss(dB, etas, dx=0.01):
    """Exact logical-error probability of a finite-energy GKP |0_D> after pure loss, decoded with the rescale rule.
    Returns dict eta -> (p_fock, sigma0_sq_measured). The unlossed state (eta=1) is simulated to measure sigma0^2."""
    Delta2 = np.arctanh(2 * 0.5 * 10 ** (-dB / 10))
    N = int(np.ceil(26 / (2 * Delta2))) + 20
    c = gkp_fock(Delta2, N); rho0 = np.outer(c, c)
    out = {}; s0 = None
    for eta in [1.0] + [e for e in etas if e != 1.0]:
        rho = rho0 if eta == 1.0 else apply_loss(rho0, eta)
        x, P = xmarginal(rho, np.sqrt(2 * N + 1) + 1, dx)
        n = np.rint(x / np.sqrt(eta) / SQRTPI)
        p = (P[(n % 2) != 0].sum() * (x[1] - x[0])) / (P.sum() * (x[1] - x[0]))
        if eta == 1.0:
            sel = np.abs(x) < SQRTPI; s0 = (P[sel] * x[sel] ** 2).sum() / P[sel].sum()
        out[eta] = (p, s0)
    return out

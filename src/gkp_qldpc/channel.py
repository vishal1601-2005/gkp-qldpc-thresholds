"""GKP displacement/loss channel layer.

Conventions (hbar = 1, [Q, P] = i, vacuum variance 1/2):
  * GKP stabilizers displace by 2*sqrt(pi); logical X/Z by sqrt(pi); decoding-cell half-width sqrt(pi)/2.
  * sigma = std-dev of the Gaussian shift per quadrature seen by the decoder on the ORIGINAL lattice.
  * Squeezing in dB: dB = -10 log10(2 sigma0^2).

Loss models (see sigma2_eff):
  rescale : loss, then a classical rescale by 1/sqrt(eta) (decoder knows eta)       sigma0^2 + (1-eta)/(2 eta)
  amp     : loss, then a physical quantum-limited amplifier of gain 1/eta           sigma0^2 + (1-eta)/eta
  paper   : the formula with NO lattice rescaling, eta*sigma0^2 + (1-eta)/2 (kept only to document that it is wrong;
            verified against an exact Fock-space simulation in gkp_qldpc.fock).
"""
import numpy as np
from scipy.special import ndtr, erfc
from scipy.optimize import brentq

SQRTPI = np.sqrt(np.pi)

def s0sq_from_dB(dB):            # sigma0^2 from squeezing in dB: dB = -10 log10(2 s0^2)
    return 0.5 * 10 ** (-dB / 10)

def sigma2_eff(s0sq, eta, model="rescale"):
    """Effective Gaussian-shift variance seen by a decoder working on the ORIGINAL lattice.
    rescale : loss, then classical rescale by 1/sqrt(eta) (decoder knows eta)  -> s0^2 + (1-eta)/(2 eta)
    amp     : loss then physical quantum-limited amplifier gain 1/eta          -> s0^2 + (1-eta)/eta
    paper   : Eq.(1) of the v1 draft (NO lattice rescaling; kept only to reproduce old numbers)
    """
    if model == "rescale": return s0sq + (1 - eta) / (2 * eta)
    if model == "amp":     return s0sq + (1 - eta) / eta
    if model == "paper":   return eta * s0sq + (1 - eta) / 2
    raise ValueError(model)

def p_flip_exact(sigma, K=40):
    """P(round(x/sqrt(pi)) odd), x~N(0,sigma^2): exact periodic sum (no nearest-cell approximation)."""
    n = np.arange(-K, K + 1)
    n = n[n % 2 != 0]
    hi = (n + 0.5) * SQRTPI / sigma
    lo = (n - 0.5) * SQRTPI / sigma
    # use survival-function form for tail accuracy: Phi(hi)-Phi(lo)
    return float(np.sum(ndtr(hi) - ndtr(lo)))

def p_flip_nearest(sigma):       # the v1 paper's Eq.(5)
    return float(erfc(np.sqrt(np.pi / 8) / sigma))

def sigma_from_p(p):             # inverse of p_flip_exact (monotone)
    return brentq(lambda s: p_flip_exact(s) - p, 1e-3, 5.0, xtol=1e-14)

def eta_threshold(sigma_th, s0sq, model="rescale"):
    """Transmissivity at which sigma2_eff == sigma_th^2 (monotone in eta)."""
    f = lambda e: sigma2_eff(s0sq, e, model) - sigma_th ** 2
    if f(1 - 1e-15) > 0: return np.nan      # squeezing alone already above threshold
    return brentq(f, 1e-9, 1 - 1e-15, xtol=1e-15)

def sample_shift_residual(sigma, size, rng):
    """x ~ N(0,sigma^2).  n=round(x/sqrt(pi)); e = n mod 2 (logical flip); r = x - n*sqrt(pi) in [-sqrt(pi)/2, sqrt(pi)/2)."""
    x = rng.normal(0.0, sigma, size)
    n = np.rint(x / SQRTPI)
    r = x - n * SQRTPI
    e = (n.astype(np.int64) % 2).astype(np.uint8)
    return r, e

def posterior_flip_prob(r, sigma, K=8):
    """P(e=1 | residual r) = sum_{n odd} g(r+n sqrt(pi)) / sum_n g(r+n sqrt(pi))  (Raveendran et al. Eq. 6)."""
    n = np.arange(-K, K + 1)
    r = np.asarray(r)[..., None]
    logg = -((r + n * SQRTPI) ** 2) / (2 * sigma ** 2)
    m = logg.max(axis=-1, keepdims=True)
    g = np.exp(logg - m)
    num = g[..., n % 2 != 0].sum(-1)
    den = g.sum(-1)
    return np.clip(num / den, 1e-12, 1 - 1e-12)

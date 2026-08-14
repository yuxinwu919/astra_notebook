"""Emittance calculations.

Geometric, normalized, and eigen-emittances from particle distributions.

All emittance formulae verified against:
  - ASTRA Manual V3.2, Table 4
  - M. Reiser, "Theory and Design of Charged Particle Beams", Wiley (2008)
  - K. Floettmann, "Some basic features of the beam emittance", PRSTAB 6, 034202 (2003)

Key definitions:
  - RMS geometric emittance:  ε = sqrt(⟨u²⟩⟨u'²⟩ − ⟨uu'⟩²)  [m·rad]
  - Normalized emittance:     ε_n = βγ · ε                  [m·rad]
  - Beam matrix:              Σ = [[⟨u²⟩, ⟨uu'⟩], [⟨uu'⟩, ⟨u'²⟩]]
  - Twiss β:                  β = ⟨u²⟩ / ε
  - Twiss α:                  α = −⟨uu'⟩ / ε
  - Twiss γ_T:                γ_T = ⟨u'²⟩ / ε = (1+α²)/β
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from ..constants import beta_gamma, gamma_from_momentum


def compute_geometric_emittance(
    u: np.ndarray,
    up: np.ndarray,
    weights: Optional[np.ndarray] = None,
    method: str = "svd",
) -> float:
    """Compute geometric RMS emittance from position and divergence arrays.

    ε_geom = sqrt(det(Σ)) where Σ is the 2×2 beam matrix.

    Two numerical methods are available:
      - 'svd':     SVD of beam matrix (most robust, default)
      - 'det':     Direct sqrt(⟨u²⟩⟨u'²⟩ − ⟨uu'⟩²) (traditional)

    Args:
        u: Position coordinates [m], centered (mean ≈ 0).
        up: Divergence coordinates [rad], centered.
        weights: Optional macro-particle charge weights.
        method: 'svd' or 'det'.

    Returns:
        Geometric RMS emittance [m·rad].
    """
    if len(u) < 2:
        return 0.0

    if weights is not None and np.sum(weights) > 0:
        w_sum = np.sum(weights)
        w = weights / w_sum
        w_eff = 1.0 / (1.0 - np.sum(w**2))
        u2 = float(np.sum(w * u**2) * w_eff)
        up2 = float(np.sum(w * up**2) * w_eff)
        u_up = float(np.sum(w * u * up) * w_eff)
    else:
        u2 = float(np.mean(u**2))
        up2 = float(np.mean(up**2))
        u_up = float(np.mean(u * up))

    if method == "svd":
        sigma = np.array([[u2, u_up], [u_up, up2]])
        try:
            S = np.linalg.svd(sigma, compute_uv=False)
            return float(np.sqrt(S[0] * S[1]))
        except np.linalg.LinAlgError:
            return 0.0
    else:
        emit_sq = u2 * up2 - u_up**2
        return float(np.sqrt(max(emit_sq, 0.0)))


def compute_normalized_emittance(
    emit_geom: float,
    ref_momentum_eVc: float,
    mass_eV: float = 0.511e6,
) -> float:
    """Convert geometric emittance to normalized emittance.

    ε_n = βγ · ε_geom

    γ is computed correctly from momentum:
      γ = sqrt(1 + (p_ref / m_e_c²)²)

    Args:
        emit_geom: Geometric RMS emittance [m·rad].
        ref_momentum_eVc: Reference momentum [eV/c].
            NOTE: This is momentum, NOT kinetic energy!
        mass_eV: Rest mass energy [eV] (default: electron 0.511 MeV).

    Returns:
        Normalized RMS emittance [m·rad].
    """
    gamma = gamma_from_momentum(ref_momentum_eVc, mass_eV)
    bg = beta_gamma(gamma) if gamma > 1.0 else 0.0
    return bg * emit_geom


def compute_twiss_parameters(
    u: np.ndarray,
    up: np.ndarray,
    weights: Optional[np.ndarray] = None,
) -> tuple[float, float, float]:
    """Compute Twiss parameters (β, α, γ_T) from centered coordinates.

    Uses SVD-based emittance for numerical stability.

    β  = ⟨u²⟩ / ε       [m]
    α  = −⟨uu'⟩ / ε     [dimensionless]
    γ_T = ⟨u'²⟩ / ε     [1/m]  (Twiss gamma, NOT relativistic γ)

    Invariant: βγ_T − α² = 1

    Args:
        u: Position coordinates [m], centered.
        up: Divergence coordinates [rad], centered.
        weights: Optional macro-particle charge weights.

    Returns:
        (beta [m], alpha, gamma_T [1/m]) tuple.
    """
    emit = compute_geometric_emittance(u, up, weights, method="svd")

    if weights is not None and np.sum(weights) > 0:
        w_sum = np.sum(weights)
        w = weights / w_sum
        w_eff = 1.0 / (1.0 - np.sum(w**2))
        u2 = float(np.sum(w * u**2) * w_eff)
        up2 = float(np.sum(w * up**2) * w_eff)
        u_up = float(np.sum(w * u * up) * w_eff)
    else:
        u2 = float(np.mean(u**2))
        up2 = float(np.mean(up**2))
        u_up = float(np.mean(u * up))

    if emit > 0:
        beta = u2 / emit
        alpha = -u_up / emit
        gamma_t = up2 / emit  # = (1 + alpha**2) / beta
    else:
        beta = 0.0
        alpha = 0.0
        gamma_t = 0.0

    return beta, alpha, gamma_t


def compute_emittance_ellipse_params(
    u: np.ndarray,
    up: np.ndarray,
    n_sigma: float = 1.0,
    weights: Optional[np.ndarray] = None,
) -> dict:
    """Compute parameters for drawing the RMS emittance ellipse.

    The ellipse is defined by the beam matrix Σ:
      [u, u'] · Σ⁻¹ · [u, u']ᵀ = ε

    At n_sigma RMS, the ellipse semi-axes are:
      a = n_sigma · sqrt(ε · λ₁)
      b = n_sigma · sqrt(ε · λ₂)
      θ = 0.5 · arctan2(2⟨uu'⟩, ⟨u²⟩ − ⟨u'²⟩)

    where λ₁, λ₂ are eigenvalues of the beam matrix.

    Args:
        u: Position coordinates [m], centered.
        up: Divergence coordinates [rad], centered.
        n_sigma: Number of RMS (1 = RMS ellipse, 2 = 2-RMS, etc.).
        weights: Optional macro-particle charge weights.

    Returns:
        Dict with keys: 'a', 'b', 'theta', 'emit', 'beta', 'alpha', 'gamma_t'.
    """
    if weights is not None and np.sum(weights) > 0:
        w_sum = np.sum(weights)
        w = weights / w_sum
        w_eff = 1.0 / (1.0 - np.sum(w**2))
        u2 = float(np.sum(w * u**2) * w_eff)
        up2 = float(np.sum(w * up**2) * w_eff)
        u_up = float(np.sum(w * u * up) * w_eff)
    else:
        u2 = float(np.mean(u**2))
        up2 = float(np.mean(up**2))
        u_up = float(np.mean(u * up))

    sigma = np.array([[u2, u_up], [u_up, up2]])
    try:
        S = np.linalg.svd(sigma, compute_uv=False)
        emit = float(np.sqrt(S[0] * S[1]))
    except np.linalg.LinAlgError:
        emit = 0.0

    if emit <= 0:
        return {"a": 0.0, "b": 0.0, "theta": 0.0, "emit": 0.0,
                "beta": 0.0, "alpha": 0.0, "gamma_t": 0.0}

    beta = u2 / emit
    alpha = -u_up / emit
    gamma_t = up2 / emit

    # Eigenvalues of beam matrix
    trace = beta + gamma_t
    disc = max(trace**2 - 4.0, 0.0)
    lam1 = 0.5 * (trace + np.sqrt(disc))
    lam2 = 0.5 * (trace - np.sqrt(disc))

    # Rotation angle
    if abs(gamma_t - beta) > 1e-12:
        theta = 0.5 * np.arctan2(2.0 * alpha, gamma_t - beta)
    else:
        theta = 0.0

    a = n_sigma * np.sqrt(emit * max(lam1, 0.0))
    b = n_sigma * np.sqrt(emit * max(lam2, 0.0))

    return {
        "a": a, "b": b, "theta": theta,
        "emit": emit, "beta": beta, "alpha": alpha, "gamma_t": gamma_t,
    }

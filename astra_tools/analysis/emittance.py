"""Emittance calculations.

Definitions follow ASTRA Manual V3.2, section 4.13, and were validated
empirically against Example.Xemit.001 (agreement better than 0.01%):

  * active particles: status > 1 (Table 2 / section 4.13)
  * canonical momentum in solenoid fields (4.13.1):
        p~_x = px + e Bz y / 2   ->   x' = p~_x / p_ref
        p~_y = py - e Bz x / 2   ->   y' = p~_y / p_ref
    where Bz is the on-axis solenoid field at the bunch center.
    In [eV/c]: p~_x = px + c Bz y / 2 (c in m/s, Bz in T, y in m).
  * geometric RMS emittance:  eps = sqrt(<u^2><u'^2> - <u u'>^2)
  * normalized emittance:     eps_n = beta*gamma * eps
    with gamma = sqrt(1 + (p_ref/mc)^2)  (momentum, NOT kinetic energy!)
  * Twiss:  beta = <u^2>/eps, alpha = -<u u'>/eps, gamma_T = <u'^2>/eps

References
----------
* ASTRA Manual V3.2, sections 4.13 and 4.13.1
* K. Floettmann, 'Some basic features of the beam emittance',
  PRST-AB 6, 034202 (2003)
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from ..constants import C_LIGHT, M_E_C2_EV, beta_gamma, gamma_from_momentum


def canonical_divergence(
    u_momentum_eVc: np.ndarray,
    v_position_m: np.ndarray,
    bz_on_axis_T: float,
    sign: float,
) -> np.ndarray:
    """Canonical transverse momentum [eV/c] for solenoid fields.

    p~_x = px + c*Bz*y/2   (sign=+1)
    p~_y = py - c*Bz*x/2   (sign=-1)

    ASTRA Manual 4.13.1: the on-axis solenoid field at the bunch center
    is used for all particles.
    """
    return u_momentum_eVc + sign * 0.5 * C_LIGHT * bz_on_axis_T * v_position_m


def compute_geometric_emittance(
    u: np.ndarray,
    up: np.ndarray,
    weights: Optional[np.ndarray] = None,
    method: str = "svd",
) -> float:
    """Geometric RMS emittance from centered coordinates [m] and
    divergence [rad].

    eps = sqrt(<u^2><u'^2> - <u u'>^2)

    Args:
        u:  position coordinates [m], centered.
        up: divergence coordinates [rad], centered.
        weights: optional macro-particle charge weights [nC]; taken in
            magnitude (|q|) so mixed-sign bunches stay weighted.
        method: 'svd' (robust, default) or 'det'.
    """
    u = np.asarray(u, dtype=float)
    up = np.asarray(up, dtype=float)
    if len(u) < 2:
        return 0.0

    if weights is not None and np.sum(np.abs(weights)) > 0:
        w = np.abs(np.asarray(weights, dtype=float))
        w = w / np.sum(w)
        # 群体矩 (ddof=0): 均匀权重与无加权分支完全一致
        u2 = float(np.sum(w * u**2))
        up2 = float(np.sum(w * up**2))
        u_up = float(np.sum(w * u * up))
    else:
        u2 = float(np.mean(u**2))
        up2 = float(np.mean(up**2))
        u_up = float(np.mean(u * up))

    if method == "svd":
        sigma = np.array([[u2, u_up], [u_up, up2]])
        try:
            s = np.linalg.svd(sigma, compute_uv=False)
            return float(np.sqrt(max(s[0] * s[1], 0.0)))
        except np.linalg.LinAlgError:
            return 0.0
    return float(np.sqrt(max(u2 * up2 - u_up**2, 0.0)))


def compute_normalized_emittance(
    emit_geom: float,
    ref_momentum_eVc: float,
    mass_eV: float = M_E_C2_EV,
) -> float:
    """Normalized emittance: eps_n = beta*gamma * eps_geom.

    gamma is computed from *momentum* [eV/c], per the ASTRA convention
    (header[1] stores p_ref).
    """
    gamma = gamma_from_momentum(ref_momentum_eVc, mass_eV)
    bg = beta_gamma(gamma)
    return bg * emit_geom


def compute_twiss_parameters(
    u: np.ndarray,
    up: np.ndarray,
    weights: Optional[np.ndarray] = None,
) -> tuple[float, float, float]:
    """Twiss parameters (beta [m], alpha, gamma_T [1/m]) from centered
    coordinates. Invariant: beta*gamma_T - alpha^2 = 1."""
    u = np.asarray(u, dtype=float)
    up = np.asarray(up, dtype=float)
    eps = compute_geometric_emittance(u, up, weights, method="svd")

    if weights is not None and np.sum(np.abs(weights)) > 0:
        w = np.abs(np.asarray(weights, dtype=float))
        w = w / np.sum(w)
        u2 = float(np.sum(w * u**2))
        up2 = float(np.sum(w * up**2))
        u_up = float(np.sum(w * u * up))
    else:
        u2 = float(np.mean(u**2))
        up2 = float(np.mean(up**2))
        u_up = float(np.mean(u * up))

    if eps > 0:
        return u2 / eps, -u_up / eps, up2 / eps
    return 0.0, 0.0, 0.0


def compute_emittance_ellipse_params(
    u: np.ndarray,
    up: np.ndarray,
    n_sigma: float = 1.0,
    weights: Optional[np.ndarray] = None,
) -> dict:
    """Parameters of the RMS emittance ellipse.

    Beam matrix Sigma = [[<u^2>, <uu'>], [<uu'>, <u'^2>]].
    Eigenvalues lambda_1 >= lambda_2 of Sigma (Twiss-normalized):
        a = n_sigma * sqrt(eps * lambda_1),  b = n_sigma * sqrt(eps * lambda_2)
        theta = 0.5 * atan2(2 alpha, gamma_T - beta)
    """
    u = np.asarray(u, dtype=float)
    up = np.asarray(up, dtype=float)
    if len(u) < 2:
        return {"a": 0.0, "b": 0.0, "theta": 0.0, "eps": 0.0,
                "beta": 0.0, "alpha": 0.0, "gamma_t": 0.0}

    if weights is not None and np.sum(np.abs(weights)) > 0:
        w = np.abs(np.asarray(weights, dtype=float))
        w = w / np.sum(w)
        u2 = float(np.sum(w * u**2))
        up2 = float(np.sum(w * up**2))
        u_up = float(np.sum(w * u * up))
    else:
        u2 = float(np.mean(u**2))
        up2 = float(np.mean(up**2))
        u_up = float(np.mean(u * up))

    sigma = np.array([[u2, u_up], [u_up, up2]])
    try:
        s = np.linalg.svd(sigma, compute_uv=False)
        eps = float(np.sqrt(max(s[0] * s[1], 0.0)))
    except np.linalg.LinAlgError:
        eps = 0.0

    if eps <= 0:
        return {"a": 0.0, "b": 0.0, "theta": 0.0, "eps": 0.0,
                "beta": 0.0, "alpha": 0.0, "gamma_t": 0.0}

    beta = u2 / eps
    alpha = -u_up / eps
    gamma_t = up2 / eps

    trace = beta + gamma_t
    disc = max(trace**2 - 4.0, 0.0)  # det = beta*gamma_t - alpha^2 = 1
    lam1 = 0.5 * (trace + np.sqrt(disc))
    lam2 = 0.5 * (trace - np.sqrt(disc))

    # 主轴角 (Floettmann): tan(2 theta) = 2<uu'> / (<u^2> - <u'^2>)
    # 不能用 Twiss 版 0.5*atan2(2*alpha, gamma_t-beta) - 它给出短轴
    # (相差 90 度), 会让 1-RMS 椭圆横躺
    if abs(u2 - up2) > 1e-30 or abs(u_up) > 1e-30:
        theta = 0.5 * np.arctan2(2.0 * u_up, u2 - up2)
    else:
        theta = 0.0

    return {
        "a": n_sigma * np.sqrt(eps * max(lam1, 0.0)),
        "b": n_sigma * np.sqrt(eps * max(lam2, 0.0)),
        "theta": theta,
        "eps": eps,
        "beta": beta,
        "alpha": alpha,
        "gamma_t": gamma_t,
    }

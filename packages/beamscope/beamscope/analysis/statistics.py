"""Beam statistical analysis.

Computes RMS sizes, emittance, Twiss parameters, energy spread, etc.
from a Distribution object.

Formulae reference:
  - RMS geometric emittance:  ε = sqrt(⟨u²⟩⟨u'²⟩ − ⟨uu'⟩²)
  - Normalized emittance:     ε_n = βγ · ε
  - Twiss β:                  β = ⟨u²⟩ / ε
  - Twiss α:                  α = −⟨uu'⟩ / ε
  - RMS (sample):             σ = sqrt(Σ(u_i − μ)² / (N−1))   [ddof=1]
  - Charge-weighted mean:     μ_w = Σ(w_i·u_i) / Σ(w_i)       [w_i = charge_i]
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional

import numpy as np

from ..constants import (
    EV_TO_MEV, M_TO_MM, EMIT_M_TO_UM,
    beta_gamma, gamma_from_momentum, kinetic_energy_from_momentum,
)
from ..distribution import Distribution


@dataclass
class BeamStatistics:
    """Comprehensive beam statistics for a particle distribution.

    All quantities are computed from active particles only.
    RMS values use sample standard deviation (ddof=1).
    """

    n_particle: int
    n_active: int
    n_lost: int = 0
    n_not_started: int = 0
    total_charge_nC: float = 0.0

    # Centroid positions [m]
    mean_x: float = 0.0
    mean_y: float = 0.0
    mean_z: float = 0.0

    # RMS sizes [m] — sample std (ddof=1)
    sig_x: float = 0.0
    sig_y: float = 0.0
    sig_z: float = 0.0

    # Mean momenta [eV/c]
    mean_px: float = 0.0
    mean_py: float = 0.0
    mean_pz: float = 0.0

    # RMS momentum spreads [eV/c] — sample std (ddof=1)
    sig_px: float = 0.0
    sig_py: float = 0.0
    sig_pz: float = 0.0

    # Relative energy spread (dimensionless)
    sig_E_over_E: float = 0.0

    # Geometric RMS emittance [m·rad]
    emit_x_geom: float = 0.0
    emit_y_geom: float = 0.0

    # Normalized RMS emittance [m·rad]
    emit_x_norm: float = 0.0
    emit_y_norm: float = 0.0

    # Twiss parameters (from RMS beam matrix)
    beta_x: float = 0.0
    alpha_x: float = 0.0
    gamma_t_x: float = 0.0
    beta_y: float = 0.0
    alpha_y: float = 0.0
    gamma_t_y: float = 0.0

    # Relativistic factors (computed from reference momentum)
    ref_momentum_eVc: float = 0.0
    ref_kinetic_energy_eV: float = 0.0
    gamma: float = 1.0
    beta_rel: float = 0.0

    # Source label
    label: str = ""

    # Charge-weighted statistics (optional, computed when use_weights=True)
    weighted: bool = False

    def to_dict(self) -> dict[str, float]:
        """Convert to a flat dictionary for pandas or serialization."""
        return {
            "n_particle": self.n_particle,
            "n_active": self.n_active,
            "n_lost": self.n_lost,
            "n_not_started": self.n_not_started,
            "total_charge_nC": self.total_charge_nC,
            "mean_x_mm": self.mean_x * M_TO_MM,
            "mean_y_mm": self.mean_y * M_TO_MM,
            "mean_z_mm": self.mean_z * M_TO_MM,
            "sig_x_mm": self.sig_x * M_TO_MM,
            "sig_y_mm": self.sig_y * M_TO_MM,
            "sig_z_mm": self.sig_z * M_TO_MM,
            "mean_pz_MeVc": self.mean_pz * EV_TO_MEV,
            "sig_E_over_E_pct": self.sig_E_over_E * 100,
            "emit_x_geom_um": self.emit_x_geom * EMIT_M_TO_UM,
            "emit_y_geom_um": self.emit_y_geom * EMIT_M_TO_UM,
            "emit_x_norm_um": self.emit_x_norm * EMIT_M_TO_UM,
            "emit_y_norm_um": self.emit_y_norm * EMIT_M_TO_UM,
            "beta_x_m": self.beta_x,
            "alpha_x": self.alpha_x,
            "gamma_t_x": self.gamma_t_x,
            "beta_y_m": self.beta_y,
            "alpha_y": self.alpha_y,
            "gamma_t_y": self.gamma_t_y,
            "ref_momentum_MeVc": self.ref_momentum_eVc * EV_TO_MEV,
            "ref_kinetic_energy_MeV": self.ref_kinetic_energy_eV * EV_TO_MEV,
            "gamma": self.gamma,
            "beta_rel": self.beta_rel,
        }


# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def compute_statistics(
    dist: Distribution,
    ref_momentum_eVc: Optional[float] = None,
    use_weights: bool = False,
    label: str = "",
) -> BeamStatistics:
    """Compute comprehensive beam statistics from a Distribution.

    All RMS quantities use sample standard deviation (ddof=1) for unbiased
    estimation from macro-particle samples.

    Normalized emittance uses γ·β computed from *momentum* (not kinetic energy),
    matching the ASTRA convention where header[1] stores p_ref [eV/c].

    Args:
        dist: Particle distribution.
        ref_momentum_eVc: Reference momentum [eV/c] for normalized emittance.
            If None, uses dist.ref_momentum_eVc or mean(pz) of active particles.
            NOTE: This is momentum, NOT kinetic energy!
        use_weights: If True, use macro-particle charge as weights for all
            statistical moments (charge-weighted mean, RMS, emittance).
        label: Optional label (e.g., 'before', 'after').

    Returns:
        BeamStatistics with all computed quantities.
    """
    mask = dist.active

    if not np.any(mask):
        raise ValueError("No active particles in distribution; cannot compute statistics.")

    x = dist.x[mask]
    y = dist.y[mask]
    z = dist.z[mask]
    px = dist.px[mask]
    py = dist.py[mask]
    pz = dist.pz[mask]
    charge = dist.charge[mask]

    # Determine reference momentum
    if ref_momentum_eVc is None:
        ref_momentum_eVc = (
            dist.ref_momentum_eVc if dist.ref_momentum_eVc != 0
            else float(np.mean(pz))
        )

    # ── Weights ──
    weights: Optional[np.ndarray] = charge if use_weights else None

    # ── Relativistic factors (from momentum, not kinetic energy!) ──
    gamma = gamma_from_momentum(ref_momentum_eVc)
    bg = beta_gamma(gamma) if gamma > 1.0 else 0.0
    beta_rel = float(np.sqrt(1.0 - 1.0 / gamma**2)) if gamma > 1.0 else 0.0
    ref_kinetic_eV = kinetic_energy_from_momentum(ref_momentum_eVc)

    # ── Helper: weighted or unweighted statistics ──
    def _mean(arr: np.ndarray, w: Optional[np.ndarray] = weights) -> float:
        if w is not None and np.sum(w) > 0:
            return float(np.average(arr, weights=w))
        return float(np.mean(arr))

    def _std(arr: np.ndarray, mu: float, w: Optional[np.ndarray] = weights) -> float:
        """Sample standard deviation (ddof=1)."""
        centered = arr - mu
        if w is not None and np.sum(w) > 0:
            w_sum = np.sum(w)
            w_norm = w / w_sum
            # Weighted sample variance: Σ w_i (x_i - μ_w)² / (1 - Σ w_i²)
            # Uses effective degrees of freedom correction
            var = np.sum(w * centered**2) / (w_sum - np.sum(w**2) / w_sum)
            return float(np.sqrt(max(var, 0.0)))
        return float(np.std(centered, ddof=1))

    # ── Centroid positions ──
    mean_x = _mean(x)
    mean_y = _mean(y)
    mean_z = _mean(z)

    # ── RMS sizes (sample std, ddof=1) ──
    sig_x = _std(x, mean_x)
    sig_y = _std(y, mean_y)
    sig_z = _std(z, mean_z)

    # ── Momentum statistics ──
    mean_px = _mean(px)
    mean_py = _mean(py)
    mean_pz = _mean(pz)

    sig_px = _std(px, mean_px)
    sig_py = _std(py, mean_py)
    sig_pz = _std(pz, mean_pz)

    sig_E_over_E = float(sig_pz / mean_pz) if mean_pz != 0 else 0.0

    # ── Transverse divergence ──
    pz_abs = np.abs(pz)
    xp = (px - mean_px) / pz_abs
    yp = (py - mean_py) / pz_abs
    x_centered = x - mean_x
    y_centered = y - mean_y

    # ── Geometric RMS emittance (SVD method, charge-weighted) ──
    emit_x_geom, beta_x, alpha_x, gamma_t_x = _compute_emittance_and_twiss(
        x_centered, xp, weights if use_weights else None,
    )
    emit_y_geom, beta_y, alpha_y, gamma_t_y = _compute_emittance_and_twiss(
        y_centered, yp, weights if use_weights else None,
    )

    # ── Normalized emittance ──
    emit_x_norm = bg * emit_x_geom if bg > 0 else 0.0
    emit_y_norm = bg * emit_y_geom if bg > 0 else 0.0

    # ── Total charge (weighted or unweighted) ──
    total_charge = float(np.sum(charge))

    return BeamStatistics(
        n_particle=dist.n_particle,
        n_active=dist.n_active,
        n_lost=int(np.sum(dist.lost)),
        n_not_started=int(np.sum(dist.not_started)),
        total_charge_nC=total_charge,
        mean_x=mean_x,
        mean_y=mean_y,
        mean_z=mean_z,
        sig_x=sig_x,
        sig_y=sig_y,
        sig_z=sig_z,
        mean_px=mean_px,
        mean_py=mean_py,
        mean_pz=mean_pz,
        sig_px=sig_px,
        sig_py=sig_py,
        sig_pz=sig_pz,
        sig_E_over_E=sig_E_over_E,
        emit_x_geom=emit_x_geom,
        emit_y_geom=emit_y_geom,
        emit_x_norm=emit_x_norm,
        emit_y_norm=emit_y_norm,
        beta_x=beta_x,
        alpha_x=alpha_x,
        gamma_t_x=gamma_t_x,
        beta_y=beta_y,
        alpha_y=alpha_y,
        gamma_t_y=gamma_t_y,
        ref_momentum_eVc=ref_momentum_eVc,
        ref_kinetic_energy_eV=ref_kinetic_eV,
        gamma=gamma,
        beta_rel=beta_rel,
        label=label,
        weighted=use_weights,
    )


def _compute_emittance_and_twiss(
    u: np.ndarray,
    up: np.ndarray,
    weights: Optional[np.ndarray] = None,
) -> tuple[float, float, float, float]:
    """Compute geometric RMS emittance and Twiss parameters via SVD.

    Uses singular value decomposition of the beam matrix for numerical
    stability, avoiding the eigenvalue discriminant issue in the
    traditional β/α/γ formulation.

    ε = sqrt(det(Σ)) where Σ = [[⟨u²⟩, ⟨uu'⟩], [⟨uu'⟩, ⟨u'²⟩]]
    β = ⟨u²⟩ / ε
    α = −⟨uu'⟩ / ε
    γ_T = ⟨u'²⟩ / ε   (Twiss gamma, distinct from relativistic γ)

    Args:
        u: Position coordinate [m], centered (mean ≈ 0).
        up: Divergence coordinate [rad], centered.
        weights: Optional macro-particle charge weights [nC].

    Returns:
        (emit_geom [m·rad], beta [m], alpha, gamma_T [1/m]) tuple.
    """
    if weights is not None and np.sum(weights) > 0:
        w_sum = np.sum(weights)
        w = weights / w_sum
        # Weighted second moments with effective DoF correction
        w_eff = 1.0 / (1.0 - np.sum(w**2))  # Kish effective sample size correction
        u2 = float(np.sum(w * u**2) * w_eff)
        up2 = float(np.sum(w * up**2) * w_eff)
        u_up = float(np.sum(w * u * up) * w_eff)
    else:
        u2 = float(np.mean(u**2))
        up2 = float(np.mean(up**2))
        u_up = float(np.mean(u * up))

    # Build 2×2 beam matrix and use SVD for robust emittance
    sigma = np.array([[u2, u_up], [u_up, up2]])
    try:
        S = np.linalg.svd(sigma, compute_uv=False)
        emit_geom = float(np.sqrt(S[0] * S[1]))
    except np.linalg.LinAlgError:
        emit_geom = 0.0

    if emit_geom > 0:
        beta = u2 / emit_geom
        alpha = -u_up / emit_geom
        gamma_t = up2 / emit_geom  # Twiss γ = (1+α²)/β
    else:
        beta = 0.0
        alpha = 0.0
        gamma_t = 0.0

    return emit_geom, beta, alpha, gamma_t


# ---------------------------------------------------------------------------
# Pretty-printer
# ---------------------------------------------------------------------------

def print_statistics(stats: BeamStatistics, title: str = "Beam Statistics") -> None:
    """Pretty-print beam statistics to stdout.

    Args:
        stats: BeamStatistics instance.
        title: Title for the output block.
    """
    width = 64
    print(f"\n{'=' * width}")
    print(f"  {title}")
    if stats.label:
        print(f"  [{stats.label}]")
    print(f"{'=' * width}")
    print(f"  Particles:       {stats.n_particle:>8d} total, "
          f"{stats.n_active:>8d} active, "
          f"{stats.n_lost:>6d} lost, "
          f"{stats.n_not_started:>6d} not started")
    print(f"  Total charge:    {stats.total_charge_nC:>10.4f} nC")
    if stats.weighted:
        print(f"  (charge-weighted statistics enabled)")
    print(f"{'-' * 40}")
    print(f"  Ref. momentum:   {stats.ref_momentum_eVc * EV_TO_MEV:>10.4f} MeV/c")
    print(f"  Ref. kinetic E:  {stats.ref_kinetic_energy_eV * EV_TO_MEV:>10.4f} MeV")
    print(f"  γ (gamma):       {stats.gamma:>10.4f}")
    print(f"  β (beta):        {stats.beta_rel:>10.6f}")
    print(f"{'-' * 40}")
    print(f"  Centroid x:      {stats.mean_x * M_TO_MM:>10.4f} mm")
    print(f"  Centroid y:      {stats.mean_y * M_TO_MM:>10.4f} mm")
    print(f"  Centroid z:      {stats.mean_z * M_TO_MM:>10.4f} mm")
    print(f"  σ_x:             {stats.sig_x * M_TO_MM:>10.4f} mm")
    print(f"  σ_y:             {stats.sig_y * M_TO_MM:>10.4f} mm")
    print(f"  σ_z:             {stats.sig_z * M_TO_MM:>10.4f} mm")
    print(f"{'-' * 40}")
    print(f"  Mean pz:         {stats.mean_pz * EV_TO_MEV:>10.4f} MeV/c")
    print(f"  σ_pz:            {stats.sig_pz * EV_TO_MEV:>10.4f} MeV/c")
    print(f"  σ_E/E:           {stats.sig_E_over_E:>10.6f}")
    print(f"{'-' * 40}")
    print(f"  ε_x (geom):      {stats.emit_x_geom * EMIT_M_TO_UM:>10.4f} μm·rad")
    print(f"  ε_y (geom):      {stats.emit_y_geom * EMIT_M_TO_UM:>10.4f} μm·rad")
    print(f"  ε_x (norm):      {stats.emit_x_norm * EMIT_M_TO_UM:>10.4f} μm·rad")
    print(f"  ε_y (norm):      {stats.emit_y_norm * EMIT_M_TO_UM:>10.4f} μm·rad")
    print(f"{'-' * 40}")
    print(f"  β_x / α_x / γ_T: {stats.beta_x:>8.4f} m / "
          f"{stats.alpha_x:>8.4f} / {stats.gamma_t_x:>8.4f} 1/m")
    print(f"  β_y / α_y / γ_T: {stats.beta_y:>8.4f} m / "
          f"{stats.alpha_y:>8.4f} / {stats.gamma_t_y:>8.4f} 1/m")
    print(f"{'=' * width}\n")

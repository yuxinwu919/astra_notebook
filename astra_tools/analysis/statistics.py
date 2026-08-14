"""Beam statistical analysis.

Computes RMS sizes, divergences, emittances, Twiss parameters and energy
spread from a Distribution. All conventions follow the ASTRA Manual V3.2
and were validated against Example.Xemit.001 (agreement < 0.01%):

  * active particles: status > 1 (manual 4.13)
  * divergence: x' = p~_x / p_ref with the canonical momentum
        p~_x = px + c Bz y / 2   (solenoid on-axis field at bunch center)
  * RMS quantities: population moments (ddof=0), matching ASTRA
  * normalized emittance: eps_n = beta*gamma * eps_geom,
    gamma = sqrt(1 + (p_ref/mc)^2)
  * energy spread: computed from per-particle kinetic energies
        E_kin = sqrt(pz^2 + m^2 c^4) - m c^2
    sigma_E/E = std(E_kin) / mean(E_kin)   (NOT sigma_p/p!)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from ..constants import (
    EMIT_M_TO_MM_MRAD,
    EV_TO_KEV,
    EV_TO_MEV,
    M_TO_MM,
    kinetic_energy_from_momentum,
    gamma_from_momentum,
    beta_from_gamma,
)
from ..distribution import Distribution
from .emittance import (
    canonical_divergence,
    compute_geometric_emittance,
    compute_twiss_parameters,
)


@dataclass
class BeamStatistics:
    """Comprehensive beam statistics (active particles only).

    SI units internally; display units in to_dict().
    """

    n_particle: int
    n_active: int
    n_passive: int = 0
    n_lost: int = 0
    n_not_started: int = 0
    total_charge_nC: float = 0.0

    # Centroids [m]
    mean_x: float = 0.0
    mean_y: float = 0.0
    mean_z: float = 0.0

    # RMS sizes [m] (population, ddof=0 - matches ASTRA)
    sig_x: float = 0.0
    sig_y: float = 0.0
    sig_z: float = 0.0

    # Mean momenta [eV/c]
    mean_px: float = 0.0
    mean_py: float = 0.0
    mean_pz: float = 0.0

    # RMS momentum spreads [eV/c]
    sig_px: float = 0.0
    sig_py: float = 0.0
    sig_pz: float = 0.0

    # Relative spreads
    sig_p_over_p: float = 0.0     # sigma_p / p
    sig_E_over_E: float = 0.0     # sigma_E / E  (kinetic energies)
    sig_E_eV: float = 0.0         # RMS kinetic-energy spread [eV]
    mean_E_kin_eV: float = 0.0    # mean kinetic energy [eV]

    # RMS divergences [rad] (canonical momentum, divided by p_ref)
    sig_xp: float = 0.0
    sig_yp: float = 0.0

    # Geometric / normalized RMS emittance [m.rad]
    emit_x_geom: float = 0.0
    emit_y_geom: float = 0.0
    emit_x_norm: float = 0.0
    emit_y_norm: float = 0.0

    # Twiss parameters
    beta_x: float = 0.0
    alpha_x: float = 0.0
    gamma_t_x: float = 0.0
    beta_y: float = 0.0
    alpha_y: float = 0.0
    gamma_t_y: float = 0.0

    # Relativistic factors from the reference momentum
    ref_momentum_eVc: float = 0.0
    ref_kinetic_energy_eV: float = 0.0
    gamma: float = 1.0
    beta_rel: float = 0.0

    # Optional canonical-momentum solenoid field
    bz_on_axis_T: float = 0.0

    label: str = ""
    weighted: bool = False

    def to_dict(self) -> dict:
        """Flat dict with accelerator-standard display units."""
        return {
            "n_particle": self.n_particle,
            "n_active": self.n_active,
            "n_passive": self.n_passive,
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
            "sig_p_over_p_pct": self.sig_p_over_p * 100,
            "mean_E_kin_MeV": self.mean_E_kin_eV * EV_TO_MEV,
            "sig_E_keV": self.sig_E_eV * EV_TO_KEV,
            "sig_E_over_E_pct": self.sig_E_over_E * 100,
            "sig_xp_mrad": self.sig_xp * 1e3,
            "sig_yp_mrad": self.sig_yp * 1e3,
            "emit_x_geom_um": self.emit_x_geom * 1e6,
            "emit_y_geom_um": self.emit_y_geom * 1e6,
            "emit_x_norm_um": self.emit_x_norm * 1e6,
            "emit_y_norm_um": self.emit_y_norm * 1e6,
            # ASTRA prints 'pi mm mrad'; numerically identical to mm.mrad
            "emit_x_norm_mm_mrad": self.emit_x_norm * EMIT_M_TO_MM_MRAD,
            "emit_y_norm_mm_mrad": self.emit_y_norm * EMIT_M_TO_MM_MRAD,
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


def compute_statistics(
    dist: Distribution,
    ref_momentum_eVc: Optional[float] = None,
    bz_on_axis_T: float = 0.0,
    use_weights: bool = False,
    label: str = "",
) -> BeamStatistics:
    """Compute comprehensive beam statistics.

    Args:
        dist: particle distribution.
        ref_momentum_eVc: reference momentum [eV/c] used for divergences
            and normalized emittance. Default: dist.ref_momentum_eVc;
            fallback to mean pz of active particles.
        bz_on_axis_T: on-axis solenoid field at the bunch center [T] for
            the canonical-momentum emittance (manual 4.13.1). Set to the
            value of the solenoid field at the bunch position if the
            emittance should be comparable to ASTRA's Xemit output.
        use_weights: charge-weighted moments (|q| weights, Kish effective
            sample size). The charge sign never silently disables the
            weighting; it only enters total_charge_nC.
        label: optional label.

    Returns:
        BeamStatistics.
    """
    mask = dist.active
    if not np.any(mask):
        raise ValueError("No active particles (status > 1) in distribution.")

    x = dist.x[mask]
    y = dist.y[mask]
    z = dist.z[mask]
    px = dist.px[mask]
    py = dist.py[mask]
    pz = dist.pz[mask]
    charge = dist.charge[mask]

    if ref_momentum_eVc is None:
        ref_momentum_eVc = (
            dist.ref_momentum_eVc if dist.ref_momentum_eVc != 0
            else float(np.mean(pz))
        )

    weights = np.abs(charge) if use_weights else None

    def _mean(a):
        if weights is not None and np.sum(weights) > 0:
            return float(np.average(a, weights=weights))
        return float(np.mean(a))

    def _std(a, mu):
        c = a - mu
        if weights is not None and np.sum(weights) > 0:
            w_sum = float(np.sum(weights))
            var = float(np.sum(weights * c**2) / (w_sum - np.sum(weights**2) / w_sum))
            return float(np.sqrt(max(var, 0.0)))
        return float(np.std(c))

    # -- Centroids / RMS sizes (population moments, matching ASTRA) --
    mean_x = _mean(x)
    mean_y = _mean(y)
    mean_z = _mean(z)
    sig_x = _std(x, mean_x)
    sig_y = _std(y, mean_y)
    sig_z = _std(z, mean_z)

    # -- Momenta --
    mean_px = _mean(px)
    mean_py = _mean(py)
    mean_pz = _mean(pz)
    sig_px = _std(px, mean_px)
    sig_py = _std(py, mean_py)
    sig_pz = _std(pz, mean_pz)
    sig_p_over_p = sig_pz / mean_pz if mean_pz != 0 else 0.0

    # -- Energy (from momentum, relativistic) --
    e_kin = kinetic_energy_from_momentum(pz)
    mean_E_kin = _mean(e_kin)
    sig_E = _std(e_kin, mean_E_kin)
    sig_E_over_E = sig_E / mean_E_kin if mean_E_kin != 0 else 0.0

    # -- Canonical divergences (manual 4.13.1) --
    ptx = canonical_divergence(px, y, bz_on_axis_T, sign=+1.0)
    pty = canonical_divergence(py, x, bz_on_axis_T, sign=-1.0)
    xp = (ptx - _mean(ptx)) / ref_momentum_eVc
    yp = (pty - _mean(pty)) / ref_momentum_eVc
    sig_xp = _std(xp, 0.0)
    sig_yp = _std(yp, 0.0)

    # -- Emittances and Twiss --
    xc = x - mean_x
    yc = y - mean_y
    emit_x_geom = compute_geometric_emittance(xc, xp, weights)
    emit_y_geom = compute_geometric_emittance(yc, yp, weights)
    beta_x, alpha_x, gamma_t_x = compute_twiss_parameters(xc, xp, weights)
    beta_y, alpha_y, gamma_t_y = compute_twiss_parameters(yc, yp, weights)

    gamma = gamma_from_momentum(ref_momentum_eVc)
    beta_rel = beta_from_gamma(gamma)
    bg = np.sqrt(max(gamma**2 - 1.0, 0.0))
    emit_x_norm = bg * emit_x_geom
    emit_y_norm = bg * emit_y_geom

    return BeamStatistics(
        n_particle=dist.n_particle,
        n_active=dist.n_active,
        n_passive=int(np.sum(dist.passive)),
        n_lost=int(np.sum(dist.lost)),
        n_not_started=int(np.sum(dist.not_started)),
        total_charge_nC=float(np.sum(charge)),
        mean_x=mean_x, mean_y=mean_y, mean_z=mean_z,
        sig_x=sig_x, sig_y=sig_y, sig_z=sig_z,
        mean_px=mean_px, mean_py=mean_py, mean_pz=mean_pz,
        sig_px=sig_px, sig_py=sig_py, sig_pz=sig_pz,
        sig_p_over_p=sig_p_over_p,
        sig_E_over_E=sig_E_over_E,
        sig_E_eV=sig_E,
        mean_E_kin_eV=mean_E_kin,
        sig_xp=sig_xp, sig_yp=sig_yp,
        emit_x_geom=emit_x_geom, emit_y_geom=emit_y_geom,
        emit_x_norm=emit_x_norm, emit_y_norm=emit_y_norm,
        beta_x=beta_x, alpha_x=alpha_x, gamma_t_x=gamma_t_x,
        beta_y=beta_y, alpha_y=alpha_y, gamma_t_y=gamma_t_y,
        ref_momentum_eVc=ref_momentum_eVc,
        ref_kinetic_energy_eV=kinetic_energy_from_momentum(ref_momentum_eVc),
        gamma=gamma, beta_rel=beta_rel,
        bz_on_axis_T=bz_on_axis_T,
        label=label,
        weighted=use_weights,
    )


def print_statistics(stats: BeamStatistics, title: str = "Beam Statistics") -> None:
    """Pretty-print beam statistics to stdout."""
    w = 64
    print("\n" + "=" * w)
    print("  " + title)
    if stats.label:
        print("  [" + stats.label + "]")
    print("=" * w)
    print("  Particles:       %8d total, %8d active, %6d passive," % (
        stats.n_particle, stats.n_active, stats.n_passive))
    print("                  %8d lost, %6d not started" % (
        stats.n_lost, stats.n_not_started))
    q_display = abs(stats.total_charge_nC)
    q_note = " (|Q| shown; sign kept internally)" if stats.total_charge_nC < 0 else ""
    print("  Total charge:    %10.4f nC%s" % (q_display, q_note))
    if stats.weighted:
        print("  (charge-weighted statistics enabled)")
    print("-" * 40)
    print("  Ref. momentum:   %10.4f MeV/c" % (stats.ref_momentum_eVc * EV_TO_MEV))
    print("  Ref. kinetic E:  %10.4f MeV" % (stats.ref_kinetic_energy_eV * EV_TO_MEV))
    print("  gamma:           %10.4f" % stats.gamma)
    print("  beta:            %10.6f" % stats.beta_rel)
    print("-" * 40)
    print("  Centroid x:      %10.4f mm" % (stats.mean_x * M_TO_MM))
    print("  Centroid y:      %10.4f mm" % (stats.mean_y * M_TO_MM))
    print("  Centroid z:      %10.4f mm" % (stats.mean_z * M_TO_MM))
    print("  sigma_x:         %10.4f mm" % (stats.sig_x * M_TO_MM))
    print("  sigma_y:         %10.4f mm" % (stats.sig_y * M_TO_MM))
    print("  sigma_z:         %10.4f mm" % (stats.sig_z * M_TO_MM))
    print("-" * 40)
    print("  Mean pz:         %10.4f MeV/c" % (stats.mean_pz * EV_TO_MEV))
    print("  sigma_p/p:       %10.6f" % stats.sig_p_over_p)
    print("  Mean E_kin:      %10.4f MeV" % (stats.mean_E_kin_eV * EV_TO_MEV))
    print("  sigma_E:         %10.4f keV" % (stats.sig_E_eV * EV_TO_KEV))
    print("  sigma_E/E:       %10.6f" % stats.sig_E_over_E)
    print("-" * 40)
    print("  eps_x (geom):    %10.4f um.rad" % (stats.emit_x_geom * 1e6))
    print("  eps_y (geom):    %10.4f um.rad" % (stats.emit_y_geom * 1e6))
    print("  eps_nx:          %10.4f um.rad (= %.4f mm.mrad)" % (
        stats.emit_x_norm * 1e6, stats.emit_x_norm * EMIT_M_TO_MM_MRAD))
    print("  eps_ny:          %10.4f um.rad (= %.4f mm.mrad)" % (
        stats.emit_y_norm * 1e6, stats.emit_y_norm * EMIT_M_TO_MM_MRAD))
    print("-" * 40)
    print("  beta_x / alpha_x / gamma_T: %8.4f m / %8.4f / %8.4f 1/m" % (
        stats.beta_x, stats.alpha_x, stats.gamma_t_x))
    print("  beta_y / alpha_y / gamma_T: %8.4f m / %8.4f / %8.4f 1/m" % (
        stats.beta_y, stats.alpha_y, stats.gamma_t_y))
    if stats.bz_on_axis_T:
        print("  (canonical momentum with Bz = %.4f T at bunch center)" % stats.bz_on_axis_T)
    print("=" * w + "\n")

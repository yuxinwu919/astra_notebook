"""Longitudinal slice analysis.

Computes slice-by-slice beam parameters (current, charge, emittance,
energy chirp) for studying slice emittance and current profiles.

Binning strategies follow the ASTRA Manual V3.2, section 6.8:
  * 'equi_spaced': equal z-width bins
  * 'equi_charge' : equal charge per bin

Physics notes (audited):
  * slices are taken from active particles (status > 1) only
  * slice divergence uses the canonical momentum (manual 4.13.1)
    normalized by the global reference momentum, u' = p~u / p_ref, and
    the normalized emittance is beta*gamma(p_ref) * eps_geom - the
    same convention as the Xemit/Yemit files
  * current: I = Q/dt with dt = dz/(beta*c); the per-slice beta/gamma
    used for the current come from the slice mean momentum
  * per-slice energy: E_kin = sqrt(pz^2 + m^2) - m (per particle)
  * equi_charge binning uses |q| (mixed-sign safe); slice charge stays
    signed internally (sign is conventional, |Q| at display)
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Optional

import numpy as np

from ..constants import C_LIGHT, NC_TO_C, kinetic_energy_from_momentum, \
    gamma_from_momentum, beta_from_gamma
from ..distribution import Distribution
from .emittance import compute_geometric_emittance


@dataclass
class SliceAnalysis:
    """Results of longitudinal slice analysis (SI units)."""

    n_slices: int
    binning_mode: str

    z_centers: np.ndarray     # [m]
    z_edges: np.ndarray       # [m] (n_slices+1)
    current: np.ndarray       # [A]
    charge: np.ndarray        # [nC]
    n_particles: np.ndarray   # int

    mean_x: np.ndarray        # [m]
    mean_y: np.ndarray        # [m]
    mean_pz: np.ndarray       # [eV/c]
    mean_kinetic_energy_eV: np.ndarray  # [eV]

    sig_x: np.ndarray         # [m]
    sig_y: np.ndarray         # [m]
    sig_E_over_E: np.ndarray  # relative energy spread

    emit_x_norm: np.ndarray   # [m.rad]
    emit_y_norm: np.ndarray   # [m.rad]

    gamma_per_slice: np.ndarray
    beta_per_slice: np.ndarray


def compute_slice_analysis(
    dist: Distribution,
    n_slices: int = 20,
    binning: str = "equi_spaced",
    ref_momentum_eVc: Optional[float] = None,
    bz_on_axis_T: float = 0.0,
) -> SliceAnalysis:
    """Compute slice-by-slice beam parameters.

    Args:
        dist: particle distribution (active particles used).
        n_slices: number of longitudinal slices.
        binning: 'equi_spaced' or 'equi_charge'.
        ref_momentum_eVc: global reference momentum [eV/c] used for the
            slice divergences and normalized emittances (Xemit convention).
        bz_on_axis_T: on-axis solenoid field at the bunch center [T] for
            the canonical momentum (manual 4.13.1).
    """
    if binning not in ("equi_spaced", "equi_charge"):
        raise ValueError("binning must be 'equi_spaced' or 'equi_charge'")

    d = dist.filter_active()
    z = d.z
    n = d.n_active
    q = d.charge

    if n < 6:
        raise ValueError("too few active particles for slice analysis")
    if n < n_slices:
        n_slices = max(n // 3, 1)

    if ref_momentum_eVc is None:
        ref_momentum_eVc = (
            d.ref_momentum_eVc if d.ref_momentum_eVc != 0
            else float(np.mean(d.pz))
        )
    if ref_momentum_eVc <= 0:
        raise ValueError(
            "reference momentum is zero/negative (beam at rest?); "
            "slice divergences and normalized emittances are undefined")

    # -- Binning --
    if binning == "equi_charge":
        sort_idx = np.argsort(z)
        z_sorted = z[sort_idx]
        q_sorted = q[sort_idx]
        q_cumsum = np.cumsum(np.abs(q_sorted))  # |q|: sign cannot break bins
        q_total = float(q_cumsum[-1])
        if q_total == 0:
            # degenerate: fall back to equi-spaced
            z_edges = np.linspace(float(np.min(z)), float(np.max(z)), n_slices + 1)
        else:
            q_per_slice = q_total / n_slices
            z_edges = np.empty(n_slices + 1)
            z_edges[0] = z_sorted[0]
            for i in range(1, n_slices):
                idx = min(np.searchsorted(q_cumsum, i * q_per_slice), n - 1)
                z_edges[i] = float(z_sorted[idx])
            z_edges[-1] = z_sorted[-1]
            # 重复 z 值会使边界塌缩 (尤其最后一条边); 用"平均箱宽
            # 的 1%" 修复全部 n_slices+1 条边并保持严格递增。delta 函数
            # 式电荷的真实峰值电流为无穷大, 这里用箱宽正则化给出有界
            # 数值 (避免 1e-12 m 级别的假 dz 造出 ~1e12 A 的假电流)
            span = float(z_sorted[-1] - z_sorted[0])
            eps = max(span / (100.0 * n_slices), 1e-12)
            for i in range(1, n_slices + 1):
                if z_edges[i] <= z_edges[i - 1]:
                    z_edges[i] = z_edges[i - 1] + eps
    else:
        z_edges = np.linspace(float(np.min(z)), float(np.max(z)), n_slices + 1)

    z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    dz = np.diff(z_edges)

    # -- Pre-allocate --
    n_part = np.zeros(n_slices, dtype=int)
    charge_s = np.zeros(n_slices)
    mean_x = np.zeros(n_slices)
    mean_y = np.zeros(n_slices)
    sig_x = np.zeros(n_slices)
    sig_y = np.zeros(n_slices)
    mean_pz = np.zeros(n_slices)
    mean_e = np.zeros(n_slices)
    sig_ee = np.zeros(n_slices)
    emit_xn = np.zeros(n_slices)
    emit_yn = np.zeros(n_slices)
    gamma_s = np.ones(n_slices)
    beta_s = np.zeros(n_slices)

    for i in range(n_slices):
        if i == n_slices - 1:
            mask = (z >= z_edges[i]) & (z <= z_edges[i + 1])
        else:
            mask = (z >= z_edges[i]) & (z < z_edges[i + 1])
        idx = np.where(mask)[0]
        n_part[i] = len(idx)
        if n_part[i] < 3:
            continue

        xi, yi = d.x[idx], d.y[idx]
        pxi, pyi, pzi = d.px[idx], d.py[idx], d.pz[idx]
        qi = q[idx]

        charge_s[i] = float(np.sum(qi))   # signed (internal); display uses |Q|
        mean_x[i] = float(np.mean(xi))
        mean_y[i] = float(np.mean(yi))
        sig_x[i] = float(np.std(xi - mean_x[i]))   # ddof=0, matches ASTRA
        sig_y[i] = float(np.std(yi - mean_y[i]))
        mean_pz[i] = float(np.mean(pzi))
        e_i = kinetic_energy_from_momentum(pzi)
        mean_e[i] = float(np.mean(e_i))
        sig_ee[i] = float(np.std(e_i) / mean_e[i]) if mean_e[i] else 0.0

        # Canonical momentum (manual 4.13.1), centered per slice, divided
        # by the global reference momentum - the Xemit convention.
        ptx = pxi + 0.5 * C_LIGHT * bz_on_axis_T * yi
        pty = pyi - 0.5 * C_LIGHT * bz_on_axis_T * xi
        xp = (ptx - np.mean(ptx)) / ref_momentum_eVc
        yp = (pty - np.mean(pty)) / ref_momentum_eVc
        # 群体矩 (无加权), 与 compute_statistics 默认及 ASTRA 一致;
        # |q| 仅用于 equi_charge 分箱
        ex = compute_geometric_emittance(xi - mean_x[i], xp)
        ey = compute_geometric_emittance(yi - mean_y[i], yp)

        # Normalized emittance w.r.t. the reference momentum (Xemit
        # convention); the per-slice gamma below is only for the current.
        g_ref = gamma_from_momentum(ref_momentum_eVc)
        bg = float(np.sqrt(max(g_ref**2 - 1.0, 0.0)))
        emit_xn[i] = bg * ex
        emit_yn[i] = bg * ey

        g_i = gamma_from_momentum(mean_pz[i])
        gamma_s[i] = g_i
        beta_s[i] = beta_from_gamma(g_i)

    # -- Current: I = Q/dt, dt = dz/(beta c) --
    current = np.zeros(n_slices)
    for i in range(n_slices):
        v = beta_s[i] * C_LIGHT if beta_s[i] > 0 else C_LIGHT
        dt = dz[i] / v if v > 0 and dz[i] > 0 else 1.0
        current[i] = charge_s[i] * NC_TO_C / dt

    return SliceAnalysis(
        n_slices=n_slices, binning_mode=binning,
        z_centers=z_centers, z_edges=z_edges,
        current=current, charge=charge_s, n_particles=n_part,
        mean_x=mean_x, mean_y=mean_y,
        mean_pz=mean_pz, mean_kinetic_energy_eV=mean_e,
        sig_x=sig_x, sig_y=sig_y, sig_E_over_E=sig_ee,
        emit_x_norm=emit_xn, emit_y_norm=emit_yn,
        gamma_per_slice=gamma_s, beta_per_slice=beta_s,
    )

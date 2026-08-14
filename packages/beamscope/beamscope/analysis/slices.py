"""Longitudinal slice analysis.

Computes slice-by-slice beam parameters for studying
energy chirp, slice emittance, and current profile.

Supports both equi-spaced and equi-charge binning strategies,
matching ASTRA's Wk_equi_grid parameter (Manual V3.2 §6.8).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal, Optional

import numpy as np

from ..constants import (
    C_LIGHT, NC_TO_C, beta_from_gamma, gamma_from_momentum,
)
from ..distribution import Distribution

BinningMode = Literal["equi_spaced", "equi_charge"]


@dataclass
class SliceAnalysis:
    """Results of longitudinal slice analysis."""

    n_slices: int
    binning_mode: BinningMode

    z_centers: np.ndarray       # [m] slice center positions
    z_edges: np.ndarray         # [m] slice boundaries (n_slices+1)
    current: np.ndarray         # [A] current per slice
    charge: np.ndarray          # [nC] charge per slice
    n_particles: np.ndarray     # particle count per slice

    # Per-slice centroids
    mean_x: np.ndarray          # [m]
    mean_y: np.ndarray          # [m]
    mean_pz: np.ndarray         # [eV/c]
    mean_kinetic_energy_eV: np.ndarray  # [eV]

    # Per-slice RMS sizes
    sig_x: np.ndarray           # [m]
    sig_y: np.ndarray           # [m]
    sig_E_over_E: np.ndarray    # relative energy spread

    # Per-slice emittance
    emit_x_norm: np.ndarray     # [m·rad] normalized x emittance
    emit_y_norm: np.ndarray     # [m·rad] normalized y emittance

    # Per-slice relativistic factors (computed from each slice's own mean pz)
    gamma_per_slice: np.ndarray
    beta_per_slice: np.ndarray


def compute_slice_analysis(
    dist: Distribution,
    n_slices: int = 20,
    binning: BinningMode = "equi_spaced",
    ref_momentum_eVc: Optional[float] = None,
) -> SliceAnalysis:
    """Compute slice-by-slice beam parameters.

    Each slice's normalized emittance uses that slice's own mean momentum
    for the βγ factor (unlike the old code which used a single reference).

    Binning strategies:
      - 'equi_spaced': Equal z-width bins (ASTRA Wk_equi_grid=1.0)
      - 'equi_charge': Equal charge per bin (ASTRA Wk_equi_grid=0.0)

    Args:
        dist: Particle distribution (only active particles used).
        n_slices: Number of longitudinal slices.
        binning: Binning strategy ('equi_spaced' or 'equi_charge').
        ref_momentum_eVc: Global reference momentum [eV/c] for fallback.
            If None, uses dist.ref_momentum_eVc.

    Returns:
        SliceAnalysis with per-slice quantities.
    """
    from .emittance import compute_geometric_emittance, compute_normalized_emittance

    d = dist.filter_active()
    z = d.z
    n = d.n_active
    charge_arr = d.charge  # already active-only

    if n < n_slices:
        n_slices = max(n // 3, 1)

    if ref_momentum_eVc is None:
        ref_momentum_eVc = (
            d.ref_momentum_eVc if d.ref_momentum_eVc != 0
            else float(np.mean(d.pz))
        )

    # ── Binning ──
    if binning == "equi_charge":
        # Sort by z, compute cumulative charge, equalize
        sort_idx = np.argsort(z)
        z_sorted = z[sort_idx]
        q_sorted = charge_arr[sort_idx]
        q_cumsum = np.cumsum(q_sorted)
        q_total = q_cumsum[-1]
        q_per_slice = q_total / n_slices

        # Find slice boundaries in sorted z
        z_edges = np.zeros(n_slices + 1)
        z_edges[0] = z_sorted[0]
        for i in range(1, n_slices):
            target_q = i * q_per_slice
            idx = np.searchsorted(q_cumsum, target_q)
            idx = min(idx, n - 1)
            z_edges[i] = float(z_sorted[idx])
        z_edges[-1] = z_sorted[-1]

        # Ensure strictly increasing (handle degenerate cases)
        for i in range(1, n_slices):
            if z_edges[i] <= z_edges[i - 1]:
                z_edges[i] = z_edges[i - 1] + 1e-12

        z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])
    else:
        # Equi-spaced
        z_min, z_max = float(np.min(z)), float(np.max(z))
        z_edges = np.linspace(z_min, z_max, n_slices + 1)
        z_centers = 0.5 * (z_edges[:-1] + z_edges[1:])

    dz_per_slice = np.diff(z_edges)

    # ── Pre-allocate ──
    n_part = np.zeros(n_slices, dtype=int)
    charge_per_slice = np.zeros(n_slices)
    mean_x = np.zeros(n_slices)
    mean_y = np.zeros(n_slices)
    sig_x = np.zeros(n_slices)
    sig_y = np.zeros(n_slices)
    mean_pz = np.zeros(n_slices)
    mean_kinetic_eV = np.zeros(n_slices)
    sig_E_over_E = np.zeros(n_slices)
    emit_x_norm = np.zeros(n_slices)
    emit_y_norm = np.zeros(n_slices)
    gamma_slice = np.ones(n_slices)
    beta_slice = np.zeros(n_slices)

    from ..constants import kinetic_energy_from_momentum

    for i in range(n_slices):
        mask = (z >= z_edges[i]) & (z < z_edges[i + 1])
        # Include last edge for the final slice
        if i == n_slices - 1:
            mask = (z >= z_edges[i]) & (z <= z_edges[i + 1])

        idx = np.where(mask)[0]
        n_part[i] = len(idx)

        if n_part[i] < 3:
            continue

        xi = d.x[idx]
        yi = d.y[idx]
        pxi = d.px[idx]
        pyi = d.py[idx]
        pzi = d.pz[idx]
        qi = charge_arr[idx]

        charge_per_slice[i] = float(np.sum(qi))
        mean_x[i] = float(np.mean(xi))
        mean_y[i] = float(np.mean(yi))
        sig_x[i] = float(np.std(xi - mean_x[i], ddof=1))
        sig_y[i] = float(np.std(yi - mean_y[i], ddof=1))
        mean_pz[i] = float(np.mean(pzi))
        mean_kinetic_eV[i] = kinetic_energy_from_momentum(mean_pz[i])

        # Per-slice divergence
        pz_abs = np.abs(pzi)
        xp = (pxi - np.mean(pxi)) / pz_abs
        yp = (pyi - np.mean(pyi)) / pz_abs

        # Per-slice emittance using this slice's own βγ
        emit_x_geom = compute_geometric_emittance(xi - mean_x[i], xp, weights=qi)
        emit_y_geom = compute_geometric_emittance(yi - mean_y[i], yp, weights=qi)

        # Each slice gets its own relativistic factor
        g_i = gamma_from_momentum(mean_pz[i])
        gamma_slice[i] = g_i
        beta_slice[i] = beta_from_gamma(g_i)
        bg_i = np.sqrt(g_i**2 - 1.0) if g_i > 1.0 else 0.0

        emit_x_norm[i] = bg_i * emit_x_geom if bg_i > 0 else 0.0
        emit_y_norm[i] = bg_i * emit_y_geom if bg_i > 0 else 0.0

        sig_pz_i = float(np.std(pzi - mean_pz[i], ddof=1))
        sig_E_over_E[i] = sig_pz_i / mean_pz[i] if mean_pz[i] != 0 else 0.0

    # ── Current: I = Q / dt, dt = dz / (βc) using per-slice β ──
    current = np.zeros(n_slices)
    for i in range(n_slices):
        v_z = beta_slice[i] * C_LIGHT if beta_slice[i] > 0 else C_LIGHT
        dt = dz_per_slice[i] / v_z if v_z > 0 and dz_per_slice[i] > 0 else 1.0
        current[i] = charge_per_slice[i] * NC_TO_C / dt

    return SliceAnalysis(
        n_slices=n_slices,
        binning_mode=binning,
        z_centers=z_centers,
        z_edges=z_edges,
        current=current,
        charge=charge_per_slice,
        n_particles=n_part,
        mean_x=mean_x,
        mean_y=mean_y,
        sig_x=sig_x,
        sig_y=sig_y,
        mean_pz=mean_pz,
        mean_kinetic_energy_eV=mean_kinetic_eV,
        sig_E_over_E=sig_E_over_E,
        emit_x_norm=emit_x_norm,
        emit_y_norm=emit_y_norm,
        gamma_per_slice=gamma_slice,
        beta_per_slice=beta_slice,
    )

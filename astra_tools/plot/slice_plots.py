"""Slice-analysis plots: current profile, slice emittance, energy chirp.

Input: astra_tools.analysis.slices.SliceAnalysis (SI); display units
applied here. Emittance display follows the ASTRA convention
([pi mm mrad] values equal eps_n in mm.mrad).
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from ..analysis.slices import SliceAnalysis


def plot_current_profile(
    sa: SliceAnalysis,
    ax=None,
    figsize=(8, 4),
    title: Optional[str] = None,
) -> plt.Figure:
    """Longitudinal current profile |I(z)| [A].

    The absolute value is shown because ASTRA stores electron macro
    charges as negative numbers (sign is conventional).
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    z = sa.z_centers * 1e3
    i_abs = np.abs(sa.current)
    ax.fill_between(z, 0, i_abs, step="mid", alpha=0.5)
    ax.plot(z, i_abs, lw=1)
    i_peak = float(np.max(i_abs))
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("current |I| [A]")
    ax.set_title(title or "longitudinal current profile")
    ax.text(0.02, 0.96, "peak %.1f A, Q = %.3f nC" % (i_peak, np.sum(sa.charge)),
            transform=ax.transAxes, va="top", fontsize=9)
    fig.tight_layout()
    return fig


def plot_slice_emittance(
    sa: SliceAnalysis,
    ax=None,
    figsize=(8, 4),
    title: Optional[str] = None,
) -> plt.Figure:
    """Slice normalized emittances [pi mm mrad] vs z."""
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    z = sa.z_centers * 1e3
    ax.plot(z, sa.emit_x_norm * 1e6, label="$\\varepsilon_{nx}$")
    ax.plot(z, sa.emit_y_norm * 1e6, label="$\\varepsilon_{ny}$")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("slice emittance [$\\pi$ mm mrad]")
    ax.set_title(title or "slice emittance")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_energy_chirp(
    sa: SliceAnalysis,
    ax=None,
    figsize=(8, 4),
    title: Optional[str] = None,
) -> plt.Figure:
    """Mean slice energy [MeV] (chirp) and relative energy spread."""
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    z = sa.z_centers * 1e3
    ax.plot(z, sa.mean_kinetic_energy_eV * 1e-6, label="$E_{kin}$", color="C0")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("mean slice $E_{kin}$ [MeV]", color="C0")
    ax2 = ax.twinx()
    ax2.plot(z, sa.sig_E_over_E * 100, label="$\\sigma_E/E$", color="C1", lw=1.2)
    ax2.set_ylabel("$\\sigma_E/E$ [%]", color="C1")
    ax.set_title(title or "energy chirp")
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [l.get_label() for l in lines], fontsize=9)
    fig.tight_layout()
    return fig


def plot_slice_sizes(
    sa: SliceAnalysis,
    ax=None,
    figsize=(8, 4),
    title: Optional[str] = None,
) -> plt.Figure:
    """Slice RMS sizes [mm] vs z."""
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    z = sa.z_centers * 1e3
    ax.plot(z, sa.sig_x * 1e3, label="$\\sigma_x$")
    ax.plot(z, sa.sig_y * 1e3, label="$\\sigma_y$")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("RMS size [mm]")
    ax.set_title(title or "slice sizes")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_slice_dashboard(
    sa: SliceAnalysis,
    figsize=(14, 8),
    title: Optional[str] = None,
) -> plt.Figure:
    """2x2 slice dashboard: current, emittance, chirp, sizes."""
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    plot_current_profile(sa, ax=axes[0, 0])
    plot_slice_emittance(sa, ax=axes[0, 1])
    plot_energy_chirp(sa, ax=axes[1, 0])
    plot_slice_sizes(sa, ax=axes[1, 1])
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig

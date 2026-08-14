"""Slice-analysis plots: current profile, slice emittance, energy chirp.

All quantities from astra_tools.analysis.slices (SI input, display units
applied here).
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
    charges as negative (the sign is conventional).
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    z = sa.z_centers * 1e3
    ax.fill_between(z, 0, np.abs(sa.current), step="mid",
                    color="steelblue", alpha=0.6)
    ax.plot(z, np.abs(sa.current), color="navy", lw=1)
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("current [A]")
    ax.set_title(title or "longitudinal current profile")
    fig.tight_layout()
    return fig


def plot_slice_emittance(
    sa: SliceAnalysis,
    ax=None,
    figsize=(8, 4),
    title: Optional[str] = None,
) -> plt.Figure:
    """Slice normalized emittances [mm.mrad] vs z."""
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    z = sa.z_centers * 1e3
    ax.plot(z, sa.emit_x_norm * 1e6, label="eps_nx", color="C0")
    ax.plot(z, sa.emit_y_norm * 1e6, label="eps_ny", color="C1")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("slice emittance [mm mrad]")
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
    ax.plot(z, sa.mean_kinetic_energy_eV * 1e-6, label="E_kin", color="C0")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("mean slice E_kin [MeV]")
    ax2 = ax.twinx()
    ax2.plot(z, sa.sig_E_over_E * 100, label="sigma_E/E", color="C1", lw=1)
    ax2.set_ylabel("sigma_E/E [%]")
    ax.set_title(title or "energy chirp")
    lines = [ax.get_lines()[0], ax2.get_lines()[0]]
    ax.legend(lines, [l.get_label() for l in lines], fontsize=8)
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

    ax = axes[1, 1]
    z = sa.z_centers * 1e3
    ax.plot(z, sa.sig_x * 1e3, label="sigma_x", color="C0")
    ax.plot(z, sa.sig_y * 1e3, label="sigma_y", color="C1")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("RMS size [mm]")
    ax.set_title("slice sizes")
    ax.legend(fontsize=8)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig

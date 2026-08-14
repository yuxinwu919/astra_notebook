"""Slice analysis visualization plots.

Slice emittance, current profile, energy chirp, and slice dashboard.
The dashboard delegates to ``PlotSliceParameters`` from the embedded
``_plotting`` core; single-panel plots are kept as lightweight wrappers.
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from ..analysis.slices import SliceAnalysis


def _slice_to_matrix(sa: SliceAnalysis) -> np.ndarray:
    """Convert SliceAnalysis to the sliceMatrix format expected by PlotSliceParameters.

    sliceMatrix columns (14):
        0=zcenter, 1=width, 2=Npart, 3=Q, 4=I,
        5=pav, 6=xcen, 7=ycen, 8=emitnx, 9=emitny,
        10=dE/E, 11=TB, 12=emitnz, 13=FullB
    """
    n = len(sa.z_centers)
    mat = np.zeros((n, 14))
    mat[:, 0] = sa.z_centers
    mat[:, 1] = np.diff(sa.z_edges) if len(sa.z_edges) == n + 1 else np.zeros(n)
    mat[:, 2] = sa.n_particles
    mat[:, 3] = sa.charge
    mat[:, 4] = sa.current
    mat[:, 5] = sa.mean_pz
    mat[:, 6] = sa.mean_x
    mat[:, 7] = sa.mean_y
    mat[:, 8] = sa.emit_x_norm
    mat[:, 9] = sa.emit_y_norm
    mat[:, 10] = sa.sig_E_over_E
    # Brightness columns (11, 13) left as 0 — not directly in SliceAnalysis
    return mat


def plot_slice_dashboard(
    sa: SliceAnalysis,
    figsize: tuple[float, float] = (14, 10),
    title: Optional[str] = None,
) -> plt.Figure:
    """Comprehensive slice analysis dashboard (3x2 grid).

    Delegates to ``PlotSliceParameters`` from the embedded ``_plotting`` core.
    """
    from .._plotting.plotter import PlotSliceParameters

    sm = _slice_to_matrix(sa)
    fig = PlotSliceParameters(sm, figsize=figsize)
    if title:
        fig.suptitle(title, fontsize=14, fontweight="bold")
    return fig


def plot_slice_emittance(
    sa: SliceAnalysis,
    figsize: tuple[float, float] = (10, 6),
    title: Optional[str] = None,
) -> plt.Figure:
    """Plot slice emittance and beam size vs z."""
    fig, ax1 = plt.subplots(1, 1, figsize=figsize)
    z_mm = sa.z_centers * 1e3

    ax1.plot(z_mm, sa.emit_x_norm * 1e6, "b-o", ms=4, lw=1.5,
             label="ε_nx [μm·rad]")
    ax1.plot(z_mm, sa.emit_y_norm * 1e6, "r--s", ms=4, lw=1.5,
             label="ε_ny [μm·rad]")
    ax1.set_xlabel("z [mm]")
    ax1.set_ylabel("Normalized Emittance [μm·rad]")
    ax1.legend(loc="upper left")
    ax1.grid(True, alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(z_mm, sa.sig_x * 1e3, "b-", lw=1.0, alpha=0.3, label="σ_x [mm]")
    ax2.plot(z_mm, sa.sig_y * 1e3, "r-", lw=1.0, alpha=0.3, label="σ_y [mm]")
    ax2.set_ylabel("RMS Size [mm]")

    ax1.set_title(title or "Slice Emittance & Beam Size")
    fig.tight_layout()
    return fig


def plot_current_profile(
    sa: SliceAnalysis,
    figsize: tuple[float, float] = (10, 4),
    title: Optional[str] = None,
    color: str = "steelblue",
) -> plt.Figure:
    """Plot longitudinal current profile."""
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    z_mm = sa.z_centers * 1e3

    ax.fill_between(z_mm, sa.current, alpha=0.6, color=color, edgecolor="white", lw=0.5)
    ax.plot(z_mm, sa.current, "k-", lw=1.0)
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("Current [A]")
    ax.set_title(title or "Longitudinal Current Profile")
    ax.grid(True, alpha=0.3)

    peak_idx = int(np.argmax(sa.current))
    peak_I = sa.current[peak_idx]
    peak_z = z_mm[peak_idx]
    ax.annotate(
        f"I_peak = {peak_I:.1f} A",
        xy=(peak_z, peak_I),
        xytext=(peak_z, peak_I * 1.1),
        ha="center", fontsize=10,
        arrowprops=dict(arrowstyle="->", color="red"),
    )

    fig.tight_layout()
    return fig


def plot_energy_chirp(
    sa: SliceAnalysis,
    figsize: tuple[float, float] = (10, 6),
    title: Optional[str] = None,
) -> plt.Figure:
    """Plot energy chirp (mean energy and energy spread per slice)."""
    fig, ax1 = plt.subplots(1, 1, figsize=figsize)
    z_mm = sa.z_centers * 1e3
    E_kin_MeV = sa.mean_kinetic_energy_eV * 1e-6

    ax1.plot(z_mm, E_kin_MeV, "b-o", ms=4, lw=1.5, label="⟨E_kin⟩ [MeV]")
    ax1.set_xlabel("z [mm]")
    ax1.set_ylabel("Mean Kinetic Energy [MeV]")
    ax1.legend(loc="upper left")
    ax1.grid(True, alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(z_mm, sa.sig_E_over_E * 100, "r--s", ms=4, lw=1.5, label="σ_E/E [%]")
    ax2.set_ylabel("Rel. Energy Spread [%]", color="r")
    ax2.tick_params(axis="y", labelcolor="r")
    ax2.legend(loc="upper right")

    ax1.set_title(title or "Energy Chirp")
    fig.tight_layout()
    return fig

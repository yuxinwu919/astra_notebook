"""Bunch Form Factor (BFF) visualization.

Delegates to ``PlotBunchFormFactor`` from the embedded ``_plotting`` core,
then adds CSR characteristic feature markers on top.
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from ..analysis.bff import BFFResult


def plot_bff(
    bff_result: BFFResult,
    x_axis: str = "wavelength",
    log_x: bool = True,
    log_y: bool = True,
    figsize: tuple[float, float] = (10, 5),
    title: Optional[str] = None,
    color: str = "steelblue",
) -> plt.Figure:
    """Plot bunch form factor |F(k)|² vs wavelength or frequency.

    Uses ``PlotBunchFormFactor`` for the base log-log plot, then adds
    CSR feature markers.
    """
    # Create figure directly (avoid wasteful create+clear from PlotBunchFormFactor)
    fig, ax = plt.subplots(figsize=figsize)
    if x_axis == "wavelength":
        x_data = bff_result.wavelength
        x_label = "Wavelength [m]"
    elif x_axis == "frequency":
        x_data = bff_result.frequency
        x_label = "Frequency [Hz]"
    else:
        x_data = bff_result.k
        x_label = "Wavenumber k [1/m]"

    ax.plot(x_data, bff_result.bff, lw=2.0, color=color, label="|F(k)|²")
    ax.set_xlabel(x_label)
    ax.set_ylabel("Bunch Form Factor |F(k)|²")

    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")

    # Mark CSR features if available
    if bff_result.csr_critical_k > 0:
        if x_axis == "wavelength":
            crit_x = 2 * np.pi / bff_result.csr_critical_k
            cutoff_x = 2 * np.pi / bff_result.csr_cutoff_k
        elif x_axis == "frequency":
            crit_x = bff_result.csr_critical_k * 2.9979e8 / (2 * np.pi)
            cutoff_x = bff_result.csr_cutoff_k * 2.9979e8 / (2 * np.pi)
        else:
            crit_x = bff_result.csr_critical_k
            cutoff_x = bff_result.csr_cutoff_k

        ax.axvline(crit_x, color="red", ls="--", lw=1.0,
                   alpha=0.7, label=f"Critical ({crit_x:.2e})")
        ax.axvline(cutoff_x, color="orange", ls=":", lw=1.0,
                   alpha=0.7, label=f"Cutoff ({cutoff_x:.2e})")

    ax.set_title(title or "Bunch Form Factor")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3, which="both")

    fig.tight_layout()
    return fig


def plot_bff_with_amplitude(
    bff_result: BFFResult,
    x_axis: str = "wavelength",
    log_x: bool = True,
    figsize: tuple[float, float] = (10, 8),
    title: Optional[str] = None,
) -> plt.Figure:
    """Plot both BFF |F(k)|² and amplitude |F(k)| in two panels."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True)

    if x_axis == "wavelength":
        x_data = bff_result.wavelength
        x_label = "Wavelength [m]"
    elif x_axis == "frequency":
        x_data = bff_result.frequency
        x_label = "Frequency [Hz]"
    else:
        x_data = bff_result.k
        x_label = "Wavenumber k [1/m]"

    # Top: BFF (squared)
    ax1.plot(x_data, bff_result.bff, lw=2.0, color="steelblue")
    ax1.set_ylabel("|F(k)|²")
    ax1.set_yscale("log")
    if log_x:
        ax1.set_xscale("log")
    ax1.grid(True, alpha=0.3, which="both")
    ax1.set_title("Bunch Form Factor (Power)")

    # Bottom: Amplitude
    ax2.plot(x_data, bff_result.bff_amplitude, lw=2.0, color="darkorange")
    ax2.set_ylabel("|F(k)|")
    ax2.set_yscale("log")
    if log_x:
        ax2.set_xscale("log")
    ax2.set_xlabel(x_label)
    ax2.grid(True, alpha=0.3, which="both")
    ax2.set_title("Bunch Form Factor (Amplitude)")
    ax2.axhline(1.0, color="gray", ls="--", lw=0.8, alpha=0.5)

    if title:
        fig.suptitle(title, fontsize=13, fontweight="bold")
    fig.tight_layout()
    return fig

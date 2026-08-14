"""Evolution plots from Xemit/Yemit/Zemit/ref files (lineplot replacement).

Standard accelerator-physics displays:
  * beam envelope: sigma_x, sigma_y [mm] vs z
  * emittance evolution: eps_nx, eps_ny [mm.mrad] vs z
  * energy: mean E_kin [MeV] and sigma_E [keV] vs z
  * reference particle: E_kin(z), dE/dz, transverse offsets, Larmor angle
  * eigen-emittances from the Sigma file (experimental, see io.astra_emit)

Units follow the ASTRA Manual V3.2 Table 4 (validated against
Example.Xemit.001).
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from ..constants import kinetic_energy_from_momentum
from ..io.astra_emit import EmitSet, RefData, SigmaData


def plot_envelope_evolution(
    emit: EmitSet,
    ax=None,
    figsize=(8, 5),
    title: Optional[str] = None,
) -> plt.Figure:
    """Beam envelope sigma_x, sigma_y [mm] vs z [m]."""
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    ax.plot(emit.x.z, emit.x.rms * 1e3, label="sigma_x", color="C0")
    ax.plot(emit.y.z, emit.y.rms * 1e3, label="sigma_y", color="C1")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("RMS beam size [mm]")
    ax.set_title(title or "beam envelope evolution")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_emittance_evolution(
    emit: EmitSet,
    ax=None,
    figsize=(8, 5),
    title: Optional[str] = None,
) -> plt.Figure:
    """Normalized emittance evolution [mm.mrad] vs z.

    Xemit/Yemit column 6 stores eps_n in units of 1e-6 m.rad
    (numerically mm.mrad) - validated against ASTRA.
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    ax.plot(emit.x.z, emit.x.emit * 1e6, label="eps_nx", color="C0")
    ax.plot(emit.y.z, emit.y.emit * 1e6, label="eps_ny", color="C1")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("normalized emittance [mm mrad]")
    ax.set_title(title or "emittance evolution")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_energy_evolution(
    emit: EmitSet,
    ax=None,
    figsize=(8, 5),
    title: Optional[str] = None,
) -> plt.Figure:
    """Mean kinetic energy [MeV] (left axis) and sigma_E [keV] (right)."""
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    ax.plot(emit.z.z, emit.z.avg * 1e-6, label="E_kin", color="C0")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("mean kinetic energy [MeV]")
    ax2 = ax.twinx()
    ax2.plot(emit.z.z, emit.z.rmsprime * 1e-3, label="sigma_E", color="C1", lw=1)
    ax2.set_ylabel("sigma_E [keV]")
    ax.set_title(title or "energy evolution")
    lines = [ax.get_lines()[0], ax2.get_lines()[0]]
    ax.legend(lines, [l.get_label() for l in lines])
    fig.tight_layout()
    return fig


def plot_ref_trajectory(
    ref: RefData,
    figsize=(12, 4),
    title: Optional[str] = None,
) -> plt.Figure:
    """Reference-particle trajectory: E_kin(z), dE/dz, offsets, Larmor."""
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    e_kin = kinetic_energy_from_momentum(ref.pz)

    ax = axes[0]
    ax.plot(ref.z, e_kin * 1e-6, color="C0")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("E_kin [MeV]")
    ax.set_title("reference energy")

    ax = axes[1]
    ax.plot(ref.z, ref.xoff * 1e3, label="x_off", color="C0")
    ax.plot(ref.z, ref.yoff * 1e3, label="y_off", color="C1")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("offset [mm]")
    ax.set_title("reference offsets")
    ax.legend(fontsize=8)

    ax = axes[2]
    ax.plot(ref.z, ref.dedz * 1e-6, color="C0")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("dE/dz [MeV/m]")
    ax.set_title("energy gain")

    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_eigen_emittances(
    sigma: SigmaData,
    ax=None,
    figsize=(8, 5),
    title: Optional[str] = None,
) -> plt.Figure:
    """Eigen-emittances from the Sigma file (EXPERIMENTAL).

    The Sigma matrix normalization is not fully documented in the manual
    (see io.astra_emit); treat as indicative only. For validated
    emittances use the Xemit files.
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    ax.plot(sigma.z, sigma.enx * 1e6, label="eps_1", color="C0")
    ax.plot(sigma.z, sigma.eny * 1e6, label="eps_2", color="C1")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("eigen-emittance [mm mrad] (experimental)")
    ax.set_title(title or "eigen-emittances (Sigma)")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_emit_dashboard(
    emit: EmitSet,
    sigma: Optional[SigmaData] = None,
    figsize=(14, 8),
    title: Optional[str] = None,
) -> plt.Figure:
    """2x2 dashboard: envelope, emittance, energy, eigen-emittances."""
    fig, axes = plt.subplots(2, 2, figsize=figsize)
    plot_envelope_evolution(emit, ax=axes[0, 0])
    plot_emittance_evolution(emit, ax=axes[0, 1])
    plot_energy_evolution(emit, ax=axes[1, 0])
    if sigma is not None:
        plot_eigen_emittances(sigma, ax=axes[1, 1])
    else:
        axes[1, 1].set_visible(False)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig

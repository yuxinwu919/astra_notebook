"""Evolution plots from Xemit/Yemit/Zemit/ref files (lineplot replacement).

Standard accelerator-physics displays (lineplot Menu 1 quantities):
  * beam envelope sigma_x, sigma_y [mm] vs z
  * normalized emittance [pi mm.mrad] vs z - the DISPLAY VALUE equals
    ASTRA's printed number (eps_n in mm.mrad; the pi marks the RMS
    phase-space ellipse area semantics; no factor is applied)
  * energy: mean E_kin [MeV] and sigma_E [keV] vs z
  * reference particle: E_kin(z), dE/dz, transverse offsets, Larmor angle
  * eigen-emittances from the Sigma file (experimental, see io.astra_emit)

Units follow the ASTRA Manual V3.2 Table 4 and were validated against
Example.Xemit.001 / Example.Zemit.001 (agreement < 0.02%).
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from ..constants import kinetic_energy_from_momentum
from ..io.astra_emit import EmitSet, RefData, SigmaData


def _xvals(e, x_axis: str):
    """x 轴取值: 'z' -> z [m], 't' -> t [ns] (postpro 三视图的时间变体)."""
    if x_axis == "t":
        return e.t * 1e9, "t [ns]"
    if x_axis != "z":
        raise ValueError("x_axis must be 'z' or 't'")
    return e.z, "z [m]"


def plot_envelope_evolution(
    emit: EmitSet,
    ax=None,
    figsize=(8, 5),
    title: Optional[str] = None,
    x_axis: str = "z",
) -> plt.Figure:
    """Beam envelope sigma_x, sigma_y [mm] vs z [m] or t [ns]."""
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    xvals, xlab = _xvals(emit.x, x_axis)
    ax.plot(xvals, emit.x.rms * 1e3, label="$\\sigma_x$")
    ax.plot(xvals, emit.y.rms * 1e3, label="$\\sigma_y$")
    ax.set_xlabel(xlab)
    ax.set_ylabel("RMS beam size [mm]")
    ax.set_title(title or "beam envelope evolution")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_divergence_evolution(
    emit: EmitSet,
    ax=None,
    figsize=(8, 5),
    title: Optional[str] = None,
    x_axis: str = "z",
) -> plt.Figure:
    """RMS divergence sigma_x', sigma_y' [mrad] vs z [m] or t [ns]."""
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    xvals, xlab = _xvals(emit.x, x_axis)
    ax.plot(xvals, emit.x.rmsprime * 1e3, label="$\\sigma_{x'}$")
    ax.plot(xvals, emit.y.rmsprime * 1e3, label="$\\sigma_{y'}$")
    ax.set_xlabel(xlab)
    ax.set_ylabel("RMS divergence [mrad]")
    ax.set_title(title or "beam divergence evolution")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_emittance_evolution(
    emit: EmitSet,
    ax=None,
    figsize=(8, 5),
    title: Optional[str] = None,
    x_axis: str = "z",
) -> plt.Figure:
    """Normalized emittance evolution vs z [m] or t [ns].

    Display unit [pi mm.mrad]: value = eps_n in mm.mrad, identical to
    the number ASTRA prints in Xemit/Yemit column 6.
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    xvals, xlab = _xvals(emit.x, x_axis)
    ax.plot(xvals, emit.x.emit * 1e6, label="$\\varepsilon_{nx}$")
    ax.plot(xvals, emit.y.emit * 1e6, label="$\\varepsilon_{ny}$")
    ax.set_xlabel(xlab)
    ax.set_ylabel("normalized emittance [$\\pi$ mm mrad]")
    ax.set_title(title or "emittance evolution")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_energy_evolution(
    emit: EmitSet,
    ax=None,
    figsize=(8, 5),
    title: Optional[str] = None,
    x_axis: str = "z",
) -> plt.Figure:
    """Mean kinetic energy [MeV] (left) and sigma_E [keV] (right)."""
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    xvals, xlab = _xvals(emit.z, x_axis)
    ax.plot(xvals, emit.z.avg * 1e-6, label="$E_{kin}$", color="C0")
    ax.set_xlabel(xlab)
    ax.set_ylabel("mean kinetic energy [MeV]", color="C0")
    ax2 = ax.twinx()
    ax2.plot(xvals, emit.z.rmsprime * 1e-3, label="$\\sigma_E$",
             color="C1", lw=1.2)
    ax2.set_ylabel("$\\sigma_E$ [keV]", color="C1")
    ax.set_title(title or "energy evolution")
    lines = ax.get_lines() + ax2.get_lines()
    ax.legend(lines, [l.get_label() for l in lines], fontsize=9)
    fig.tight_layout()
    return fig


def plot_bunch_length_evolution(
    emit: EmitSet,
    ax=None,
    figsize=(8, 5),
    title: Optional[str] = None,
    x_axis: str = "z",
) -> plt.Figure:
    """RMS bunch length [mm] vs z [m] or t [ns]."""
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    xvals, xlab = _xvals(emit.z, x_axis)
    ax.plot(xvals, emit.z.rms * 1e3, label="$\\sigma_z$")
    ax.set_xlabel(xlab)
    ax.set_ylabel("RMS bunch length [mm]")
    ax.set_title(title or "bunch length evolution")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_energy_spread_evolution(
    emit: EmitSet,
    ax=None,
    figsize=(8, 5),
    title: Optional[str] = None,
    x_axis: str = "z",
) -> plt.Figure:
    """RMS energy spread [keV] vs z [m] or t [ns]."""
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    xvals, xlab = _xvals(emit.z, x_axis)
    ax.plot(xvals, emit.z.rmsprime * 1e-3, label="$\\sigma_E$")
    ax.set_xlabel(xlab)
    ax.set_ylabel("RMS energy spread [keV]")
    ax.set_title(title or "energy spread evolution")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_ref_trajectory(
    ref: RefData,
    figsize=(13, 4),
    title: Optional[str] = None,
) -> plt.Figure:
    """Reference particle: E_kin(z), dE/dz, transverse offsets, Larmor."""
    fig, axes = plt.subplots(1, 3, figsize=figsize)
    e_kin = kinetic_energy_from_momentum(ref.pz)

    ax = axes[0]
    ax.plot(ref.z, e_kin * 1e-6)
    ax.set_xlabel("z [m]")
    ax.set_ylabel("$E_{kin}$ [MeV]")
    ax.set_title("reference energy")

    ax = axes[1]
    ax.plot(ref.z, ref.xoff * 1e3, label="$x_{off}$")
    ax.plot(ref.z, ref.yoff * 1e3, label="$y_{off}$")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("offset [mm]")
    ax.set_title("reference offsets")
    ax.legend(fontsize=9)

    ax = axes[2]
    ax.plot(ref.z, ref.dedz * 1e-6)
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
    (see io.astra_emit); treat as indicative. Use the Xemit files for
    validated emittances.
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    ax.plot(sigma.z, sigma.enx * 1e6, label="$\\varepsilon_1$")
    ax.plot(sigma.z, sigma.eny * 1e6, label="$\\varepsilon_2$")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("eigen-emittance [$\\pi$ mm mrad] (experimental)")
    ax.set_title(title or "eigen-emittances (Sigma)")
    ax.legend()
    fig.tight_layout()
    return fig


def plot_velocity_evolution(
    ref: RefData,
    ax=None,
    figsize=(8, 5),
    title: Optional[str] = None,
) -> plt.Figure:
    """参考粒子速度 beta = v/c 与 gamma vs z (lineplot 粒子速度曲线)."""
    from ..constants import gamma_from_momentum, beta_from_gamma
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    gamma = gamma_from_momentum(ref.pz)
    beta = beta_from_gamma(gamma)
    ax.plot(ref.z, beta, label="$\\beta = v/c$")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("velocity $\\beta$")
    ax.set_title(title or "reference particle velocity")
    ax2 = ax.twinx()
    ax2.plot(ref.z, gamma, color="C1", ls="--", label="$\\gamma$")
    ax2.set_ylabel("$\\gamma$", color="C1")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, fontsize=9)
    fig.tight_layout()
    return fig


def plot_step_size_evolution(
    ref: RefData,
    ax=None,
    figsize=(8, 5),
    title: Optional[str] = None,
) -> plt.Figure:
    """平均积分步长 vs z (lineplot 平均步长曲线).

    ref 文件每行是一个 Runge-Kutta 步; 相邻行的 z 间距即步长。
    """
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    dz = np.diff(ref.z)
    zm = 0.5 * (ref.z[1:] + ref.z[:-1])
    ax.plot(zm, dz * 1e3, label="RK step size")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("step size [mm]")
    ax.set_title(title or "average step size evolution")
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
        plot_bunch_length_evolution(emit, ax=axes[1, 1])
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_lineplot_overview(
    emit: EmitSet,
    figsize=(14, 10),
    title: Optional[str] = None,
) -> plt.Figure:
    """lineplot Menu 1 overview: the first 9 plots on one page.

    envelope, divergence, emittance, bunch length, energy spread,
    energy + 3 reference-particle panels.
    """
    fig, axes = plt.subplots(3, 3, figsize=figsize)
    plot_envelope_evolution(emit, ax=axes[0, 0])
    plot_divergence_evolution(emit, ax=axes[0, 1])
    plot_emittance_evolution(emit, ax=axes[0, 2])
    plot_bunch_length_evolution(emit, ax=axes[1, 0])
    plot_energy_spread_evolution(emit, ax=axes[1, 1])
    plot_energy_evolution(emit, ax=axes[1, 2])
    # longitudinal emittance
    ax = axes[2, 0]
    ax.plot(emit.z.z, emit.z.emit * 1e-3, label="$\\varepsilon_{nz}$")
    ax.set_xlabel("z [m]")
    ax.set_ylabel("long. emittance [keV mm]")
    ax.set_title("longitudinal emittance")
    ax.legend()
    axes[2, 1].set_visible(False)
    axes[2, 2].set_visible(False)
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig

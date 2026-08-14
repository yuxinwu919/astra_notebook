"""ASTRA evolution plots: envelope, emittance, energy, eigen-emittances, reference.

All functions take ``EmitSet``, ``SigmaData``, or ``RefData`` and delegate
to the embedded ``_plotting`` core for battle-tested matplotlib rendering.
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from beamscope.io.astra_emit import EmitSet, EmitData, SigmaData, RefData
from beamscope._plotting.plotter import (
    PlotSize1plt,
    PlotSize1pltLat,
    PlotTransSize1plt,
    PlotTransSize1pltMag,
    PlotEmit1plt,
    PlotEnergy1plt,
    PlotEigenEmits,
)

# ASTRA emit dtype used by the _plotting core
_EMIT_DTYPE = np.dtype({
    'names': ['z', 't', 'avg', 'rms', 'rmsprime', 'emit', 'corr'],
    'formats': [np.float64] * 7
})


def _emitdata_to_array(d: EmitData) -> np.ndarray:
    """Convert an EmitData dataclass to a structured numpy array."""
    n = len(d.z)
    arr = np.zeros(n, dtype=_EMIT_DTYPE)
    arr['z'] = d.z
    arr['t'] = d.t
    arr['avg'] = d.avg
    arr['rms'] = d.rms
    arr['rmsprime'] = d.rmsprime
    arr['emit'] = d.emit
    arr['corr'] = d.corr
    return arr


def _fig_or_new(fig: Optional[plt.Figure] = None, figsize=(10, 5)) -> plt.Figure:
    """Return fig if provided, otherwise create a new one."""
    if fig is not None:
        fig.clear()
        return fig
    return plt.figure(figsize=figsize)


# ============================================================
# Plot functions — delegate to _plotting core
# ============================================================

def plot_envelope_evolution(
    emit: EmitSet,
    fig: Optional[plt.Figure] = None,
    figsize: tuple = (10, 5),
    title: Optional[str] = None,
) -> plt.Figure:
    """Plot beam envelope (RMS sizes) vs z.

    Delegates to ``PlotSize1plt`` from the embedded plotting core.
    """
    X = _emitdata_to_array(emit.x)
    Y = _emitdata_to_array(emit.y)
    Z = _emitdata_to_array(emit.z)
    fig = PlotSize1plt(X, Y, Z, figsize=figsize)
    if title:
        fig.axes[0].set_title(title)
    elif emit.filename:
        fig.axes[0].set_title(f"Beam Envelope — {emit.filename}")
    return fig


def plot_transverse_size(
    emit: EmitSet,
    fig: Optional[plt.Figure] = None,
    figsize: tuple = (10, 5),
    title: Optional[str] = None,
) -> plt.Figure:
    """Plot transverse RMS beam sizes only (no longitudinal).

    Delegates to ``PlotTransSize1plt`` from the embedded plotting core.
    """
    X = _emitdata_to_array(emit.x)
    Y = _emitdata_to_array(emit.y)
    fig = PlotTransSize1plt(X, Y, figsize=figsize)
    if title:
        fig.axes[0].set_title(title)
    elif emit.filename:
        fig.axes[0].set_title(f"Transverse Size — {emit.filename}")
    return fig


def plot_size_with_lattice(
    emit: EmitSet,
    lattice: np.ndarray,
    fig: Optional[plt.Figure] = None,
    figsize: tuple = (12, 6),
    title: Optional[str] = None,
) -> plt.Figure:
    """Plot beam sizes with lattice profile overlay.

    Parameters
    ----------
    emit : EmitSet
    lattice : ndarray
        Structured array with 'z' and 'profile' fields for the lattice.
    """
    X = _emitdata_to_array(emit.x)
    Y = _emitdata_to_array(emit.y)
    Z = _emitdata_to_array(emit.z)
    fig = PlotSize1pltLat(X, Y, Z, lattice, figsize=figsize)
    if title:
        fig.axes[0].set_title(title)
    elif emit.filename:
        fig.axes[0].set_title(f"Size + Lattice — {emit.filename}")
    return fig


def plot_size_with_magnets(
    emit: EmitSet,
    magnets: np.ndarray,
    mag_offset: float = 0.0,
    mag_scale: float = 1.0,
    fig: Optional[plt.Figure] = None,
    figsize: tuple = (10, 6),
    title: Optional[str] = None,
) -> plt.Figure:
    """Plot transverse beam sizes with magnet profile overlay.

    Parameters
    ----------
    emit : EmitSet
    magnets : ndarray
        Structured array with 'z' and 'profile' fields.
    mag_offset : float
        Vertical offset for magnet profile.
    mag_scale : float
        Scale factor for magnet profile.
    """
    X = _emitdata_to_array(emit.x)
    Y = _emitdata_to_array(emit.y)
    fig = PlotTransSize1pltMag(X, Y, magnets, mag_offset, mag_scale, figsize=figsize)
    if title:
        fig.axes[0].set_title(title)
    elif emit.filename:
        fig.axes[0].set_title(f"Size + Magnets — {emit.filename}")
    return fig


def plot_emittance_evolution(
    emit: EmitSet,
    fig: Optional[plt.Figure] = None,
    figsize: tuple = (10, 5),
    title: Optional[str] = None,
) -> plt.Figure:
    """Plot normalized emittance evolution vs z.

    Delegates to ``PlotEmit1plt`` from the embedded plotting core.
    """
    X = _emitdata_to_array(emit.x)
    Y = _emitdata_to_array(emit.y)
    Z = _emitdata_to_array(emit.z)
    fig = PlotEmit1plt(X, Y, Z, figsize=figsize)
    if title:
        fig.axes[0].set_title(title)
    elif emit.filename:
        fig.axes[0].set_title(f"Normalized Emittance — {emit.filename}")
    return fig


def plot_energy_evolution(
    emit: EmitSet,
    fig: Optional[plt.Figure] = None,
    figsize: tuple = (10, 5),
    title: Optional[str] = None,
) -> plt.Figure:
    """Plot energy and energy spread vs z.

    Delegates to ``PlotEnergy1plt`` from the embedded plotting core.
    """
    X = _emitdata_to_array(emit.x)
    Y = _emitdata_to_array(emit.y)
    Z = _emitdata_to_array(emit.z)
    fig = PlotEnergy1plt(X, Y, Z, figsize=figsize)
    if title:
        fig.axes[0].set_title(title)
    elif emit.filename:
        fig.axes[0].set_title(f"Energy Evolution — {emit.filename}")
    return fig


def plot_eigen_emittances(
    sigma: SigmaData,
    fig: Optional[plt.Figure] = None,
    figsize: tuple = (10, 5),
    title: Optional[str] = None,
) -> plt.Figure:
    """Plot eigen-emittances from sigma matrix vs z.

    Delegates to ``PlotEigenEmits`` from the embedded plotting core.
    """
    # Build a minimal structured array with 'z' field for the core plotter
    S = np.zeros(len(sigma.z), dtype=[('z', 'f8')])
    S['z'] = sigma.z
    fig = PlotEigenEmits(S, sigma.enx, sigma.eny, sigma.enz, figsize=figsize)
    if title:
        fig.axes[0].set_title(title)
    elif sigma.filename:
        fig.axes[0].set_title(f"Eigen-Emittances — {sigma.filename}")
    return fig


def plot_ref_trajectory(
    ref: RefData,
    fig: Optional[plt.Figure] = None,
    figsize: tuple = (10, 8),
    title: Optional[str] = None,
) -> plt.Figure:
    """Plot reference particle trajectory: position, momentum, energy vs z.

    2x2 subplots. (Custom — no direct astra_plotter equivalent.)
    """
    fig = _fig_or_new(fig, figsize)
    z = ref.z

    axes = fig.subplots(2, 2)

    # Position
    ax = axes[0, 0]
    ax.plot(z, ref.xoff, "b-", lw=1.5, label="xoff")
    ax.plot(z, ref.yoff, "r--", lw=1.5, label="yoff")
    ax.set_xlabel("z [m]"); ax.set_ylabel("Offset [mm]")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    ax.set_title("Reference Position")

    # Transverse momentum
    ax = axes[0, 1]
    ax.plot(z, ref.px, "b-", lw=1.5, label="px")
    ax.plot(z, ref.py, "r--", lw=1.5, label="py")
    ax.set_xlabel("z [m]"); ax.set_ylabel("Momentum [eV/c]")
    ax.legend(fontsize=8); ax.grid(True, alpha=0.3)
    ax.set_title("Transverse Momentum")

    # Longitudinal momentum
    ax = axes[1, 0]
    ax.plot(z, ref.pz, "g-", lw=1.5)
    ax.set_xlabel("z [m]"); ax.set_ylabel("pz [MeV/c]")
    ax.grid(True, alpha=0.3)
    ax.set_title("Longitudinal Momentum")

    # Energy loss
    ax = axes[1, 1]
    ax.plot(z, ref.dE_dz, "purple", lw=1.5)
    ax.set_xlabel("z [m]"); ax.set_ylabel("dE/dz [MeV/m]")
    ax.grid(True, alpha=0.3)
    ax.set_title("Energy Gradient")

    fig.suptitle(title or f"Reference Particle — {ref.filename}", fontweight="bold")
    try:
        fig.tight_layout()
    except Exception:
        pass
    return fig


def plot_emit_dashboard(
    emit: EmitSet,
    sigma: Optional[SigmaData] = None,
    fig: Optional[plt.Figure] = None,
    figsize: tuple = (14, 10),
    title: Optional[str] = None,
) -> plt.Figure:
    """Comprehensive emit dashboard: envelope + emittance + energy + eigen-emittances.

    2x2 grid. Each panel delegates to the embedded ``_plotting`` core for
    inline rendering (no per-panel figure creation — we reuse subplots).
    """
    fig = _fig_or_new(fig, figsize)
    axes = fig.subplots(2, 2)

    # Envelope (top-left): inline rendering for subplot embedding
    ax = axes[0, 0]
    z = emit.x.z
    ax.plot(z, emit.x.rms * 1e3, "b-", lw=1.2, label="sig_x")
    ax.plot(z, emit.y.rms * 1e3, "r--", lw=1.2, label="sig_y")
    ax2_e = ax.twinx()
    ax2_e.plot(z, emit.z.rms * 1e3, "g-.", lw=1.2, label="sig_z")
    ax.set_xlabel("z [m]"); ax.set_ylabel("sig_x, sig_y [mm]")
    ax2_e.set_ylabel("sig_z [mm]", color="g")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(True, alpha=0.3); ax.set_title("Beam Envelope")

    # Emittance (top-right)
    ax = axes[0, 1]
    ax.plot(z, emit.x.emit, "b-", lw=1.2, label="g*eps_x")
    ax.plot(z, emit.y.emit, "r--", lw=1.2, label="g*eps_y")
    ax2_em = ax.twinx()
    ax2_em.plot(z, emit.z.emit, "g-.", lw=1.2, label="g*eps_z")
    ax.set_xlabel("z [m]"); ax.set_ylabel("g*eps_x, g*eps_y [pi*mrad*mm]")
    ax2_em.set_ylabel("g*eps_z [pi*keV*mm]", color="g")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(True, alpha=0.3); ax.set_title("Normalized Emittance")

    # Energy (bottom-left)
    ax = axes[1, 0]
    sig_E = emit.z.rmsprime * 1e-3
    corr_E = np.abs(emit.z.corr) * 1e-3
    E_kin = emit.z.avg * 1e-6
    ax.plot(z, sig_E, "b-", lw=1.2, label="sig_E RMS")
    ax.plot(z, corr_E, "r--", lw=1.2, label="corr sig_E")
    ax2_en = ax.twinx()
    ax2_en.plot(z, E_kin, "g-", lw=1.2)
    ax.set_xlabel("z [m]"); ax.set_ylabel("Energy Spread [keV]")
    ax2_en.set_ylabel("E_kin [MeV]", color="g")
    ax.legend(fontsize=7, loc="upper left")
    ax.grid(True, alpha=0.3); ax.set_title("Energy Evolution")

    # Eigen-emittances (bottom-right)
    ax = axes[1, 1]
    if sigma is not None:
        ax.plot(sigma.z, sigma.enx * 1e6, "b-", lw=1.5, label="eps_nx")
        ax.plot(sigma.z, sigma.eny * 1e6, "r--", lw=1.5, label="eps_ny")
        ax2_ee = ax.twinx()
        ax2_ee.plot(sigma.z, sigma.enz * 1e6, "g-.", lw=1.5, label="eps_nz")
        ax.set_xlabel("z [m]"); ax.set_ylabel("eps_nx, eps_ny [um*rad]")
        ax2_ee.set_ylabel("eps_nz [um]", color="g")
        ax2_ee.tick_params(axis="y", labelcolor="g")
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)
        ax.set_title("Eigen-Emittances")
    else:
        ax.text(0.5, 0.5, "No Sigma file loaded\n\nLoad .Sigma file for\neigen-emittance analysis",
                transform=ax.transAxes, ha="center", va="center", fontsize=12, color="gray")
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title("Eigen-Emittances (unavailable)")

    fig.suptitle(title or f"Emit Dashboard — {emit.filename}", fontsize=14, fontweight="bold")
    try:
        fig.tight_layout()
    except Exception:
        pass
    return fig

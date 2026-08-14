"""
ASTRA Plotter — High-level plotting functions for beam dynamics data.

Provides publication-quality plotting functions for:
- Emittance evolution (transverse & longitudinal)
- RMS beam size evolution
- Energy and energy spread
- Phase space density plots (2D histograms with projections)
- Eigen-emittance plots
- Bunch form factor plots
- Slice parameter dashboards

All functions return matplotlib Figure objects.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

from .cosmetics import FormatLabelSci, SLAC_DESY_CMAP

beam_map = SLAC_DESY_CMAP  # legacy alias


# ============================================================================
# Helpers
# ============================================================================

def histogram0(x, Nbins):
    """Create a histogram with zero-valued boundary bins (for step-plot rendering).

    Parameters
    ----------
    x : ndarray
        Input data.
    Nbins : int
        Number of bins.

    Returns
    -------
    yhist0 : ndarray
        Histogram counts (Nbins+2), zero-padded.
    xhist0 : ndarray
        Bin edges (Nbins+2).
    """
    yhist, xhist = np.histogram(x, bins=Nbins)
    yhist0 = np.zeros(Nbins + 2)
    xhist0 = np.zeros(Nbins + 2)
    xhist0[1:Nbins + 1] = xhist[0:Nbins]
    yhist0[1:Nbins + 1] = yhist[0:Nbins]
    dx = xhist[1] - xhist[0]
    xhist0[0] = xhist[0] - dx
    xhist0[Nbins + 1] = xhist[Nbins - 1] + dx
    return (yhist0, xhist0)


# ============================================================================
# Helper
# ============================================================================

def _add_ymargin_if_flat(ax, data, threshold=0.001, margin=0.05):
    """Add y-axis margin if data has very small relative variation."""
    dmin, dmax = np.min(data), np.max(data)
    if abs(dmax) < 1e-30:
        return
    relative_variation = (dmax - dmin) / abs(dmax)
    if relative_variation < threshold:
        center = (dmin + dmax) / 2.0
        half_span = max(abs(dmax - center), abs(center - dmin)) + abs(center) * margin
        ax.set_ylim(center - half_span, center + half_span)


# ============================================================================
# Emittance Plots
# ============================================================================

def PlotEmit1plt(X, Y, Z, figsize=(10, 6)):
    """Plot transverse and longitudinal emittance evolution on dual y-axes.

    Parameters
    ----------
    X, Y : ndarray
        Transverse emit structured arrays (fields: 'z', 'emit').
    Z : ndarray
        Longitudinal emit structured array.
    figsize : tuple
        Figure size.
    """
    fig, ax1 = plt.subplots(figsize=figsize)

    ax1.plot(X['z'], X['emit'], '-', color='blue', linewidth=2.0,
             label=r'$\gamma \epsilon_{x}$')
    ax1.plot(Y['z'], Y['emit'], '--', color='red', linewidth=2.0,
             label=r'$\gamma \epsilon_{y}$')
    ax1.legend(loc='lower right')
    ax1.set_ylabel(r'Transverse emittance $\gamma\epsilon_{x,y}$ ($\mu$m)', fontsize=22)
    ax1.set_xlabel(r'Distance $z$ (m)', fontsize=22)
    _add_ymargin_if_flat(ax1, np.concatenate([X['emit'], Y['emit']]))
    FormatLabelSci()

    ax2 = ax1.twinx()
    ax2.plot(Z['z'], Z['emit'], color='green', linewidth=2.0,
             label=r'$\gamma \epsilon_{z}$')
    ax2.set_ylabel(r'Longitudinal emittance $\gamma\epsilon_{z}$ ($\mu$m)', fontsize=22,
                   color="green")
    FormatLabelSci()
    for label in ax2.get_yticklabels():
        label.set_color("green")
    ax2.legend(loc='upper left')
    plt.tight_layout()
    return fig


def PlotEigenEmits(S, enx, eny, enz, figsize=(10, 6)):
    """Plot eigen-emittances vs longitudinal position.

    Parameters
    ----------
    S : ndarray
        Sigma structured array (field: 'z').
    enx, eny : ndarray
        Transverse eigen-emittances.
    enz : ndarray
        Longitudinal emittance.
    figsize : tuple
        Figure size.
    """
    fig, ax1 = plt.subplots(figsize=figsize)

    ax1.plot(S['z'], enx, '-', color='blue', linewidth=2.0,
             label=r'$\gamma \epsilon_{1}$')
    ax1.plot(S['z'], eny, '--', color='red', linewidth=2.0,
             label=r'$\gamma \epsilon_{2}$')
    ax1.legend(loc='lower right')
    ax1.set_ylabel(r'Transverse eigen-emittance $\gamma\epsilon_{1,2}$ ($\mu$m)', fontsize=22)
    ax1.set_xlabel(r'Distance $z$ (m)', fontsize=22)
    FormatLabelSci()

    ax2 = ax1.twinx()
    ax2.plot(S['z'], enz, color='green', linewidth=2.0,
             label=r'$\gamma \epsilon_{z}$')
    ax2.set_ylabel(r'Longitudinal emittance $\gamma\epsilon_{z}$ ($\mu$m)', fontsize=22,
                   color="green")
    FormatLabelSci()
    for label in ax2.get_yticklabels():
        label.set_color("green")
    ax2.legend(loc='upper left')
    plt.tight_layout()
    return fig


# ============================================================================
# Beam Size Plots
# ============================================================================

def PlotSize1plt(X, Y, Z, figsize=(10, 6)):
    """Plot transverse and longitudinal RMS beam sizes on dual y-axes.

    Parameters
    ----------
    X, Y : ndarray
        Transverse emit structured arrays (fields: 'z', 'rms').
    Z : ndarray
        Longitudinal emit structured array.
    figsize : tuple
        Figure size.
    """
    fig, ax1 = plt.subplots(figsize=figsize)

    ax1.plot(X['z'], X['rms'], color='blue', linewidth=2.5,
             label=r'$\sigma_{x}$ (mm)')
    ax1.plot(Y['z'], Y['rms'], '--', color='red', linewidth=2.0,
             label=r'$\sigma_{y}$ (mm)')
    ax1.set_ylabel(r'Transverse RMS beam size (mm)', fontsize=22)
    ax1.set_xlabel(r'Distance $z$ (m)', fontsize=22)
    ax1.legend(loc='upper left')
    _add_ymargin_if_flat(ax1, X['rms'])
    FormatLabelSci()

    ax2 = ax1.twinx()
    ax2.plot(Z['z'], Z['rms'], '-.', color='green', linewidth=2.0,
             label=r'$\sigma_{z}$ (mm)')
    ax2.set_ylabel(r'Longitudinal RMS bunch length (mm)', fontsize=22,
                   color="green")
    ax2.legend(loc='upper right')
    _add_ymargin_if_flat(ax2, Z['rms'])
    for label in ax2.get_yticklabels():
        label.set_color("green")
    FormatLabelSci()
    plt.tight_layout()
    return fig


def PlotSize1pltLat(X, Y, Z, Latt, figsize=(12, 6)):
    """Plot beam sizes with lattice profile overlay.

    Parameters
    ----------
    X, Y, Z : ndarray
        Emit arrays.
    Latt : ndarray
        Lattice profile with 'z' and 'profile' fields.
    figsize : tuple
        Figure size.
    """
    fig, ax1 = plt.subplots(figsize=figsize)

    ax1.plot(X['z'], X['rms'], color='blue', linewidth=2.5,
             label=r'$\sigma_{x}$ (mm)')
    ax1.plot(Y['z'], Y['rms'], '--', color='red', linewidth=2.0,
             label=r'$\sigma_{y}$ (mm)')
    ax1.plot(Latt['z'], Latt['profile'], '--', color='grey', linewidth=2.0,
             label='Lattice')
    ax1.legend(loc='upper left')
    ax1.set_ylabel(r'Transverse RMS beam size (mm)', fontsize=22)
    ax1.set_xlabel(r'Distance $z$ (m)', fontsize=22)
    FormatLabelSci()

    ax2 = ax1.twinx()
    ax2.plot(Z['z'], Z['rms'], '-.', color='green', linewidth=2.0,
             label=r'$\sigma_{z}$ (mm)')
    ax2.set_ylabel(r'Longitudinal RMS bunch length (mm)', fontsize=22,
                   color="green")
    ax2.legend(loc='upper right')
    for label in ax2.get_yticklabels():
        label.set_color("green")
    FormatLabelSci()
    plt.tight_layout()
    return fig


def PlotTransSize1plt(X, Y, figsize=(10, 5)):
    """Plot transverse RMS beam sizes only.

    Parameters
    ----------
    X, Y : ndarray
        Transverse emit arrays (fields: 'z', 'rms').
    figsize : tuple
        Figure size.
    """
    fig, ax1 = plt.subplots(figsize=figsize)

    ax1.plot(X['z'], X['rms'], color='blue', linewidth=2.5,
             label=r'$\sigma_{x}$')
    ax1.plot(Y['z'], Y['rms'], '--', color='red', linewidth=2.0,
             label=r'$\sigma_{y}$')
    ax1.legend()
    ax1.set_ylabel(r'Transverse RMS beam size (mm)', fontsize=22)
    ax1.set_xlabel(r'Distance $z$ (m)', fontsize=22)
    _add_ymargin_if_flat(ax1, np.concatenate([X['rms'], Y['rms']]))
    FormatLabelSci()
    plt.tight_layout()
    return fig


def PlotTransSize1pltMag(X, Y, MAG, MAGoffset=0., MAGscale=1.,
                         figsize=(10, 6)):
    """Plot transverse beam sizes with magnet profile overlay.

    Parameters
    ----------
    X, Y : ndarray
        Transverse emit arrays.
    MAG : ndarray
        Magnet profile with 'z' and 'profile' fields.
    MAGoffset, MAGscale : float
        Offset and scale for magnet profile.
    figsize : tuple
        Figure size.
    """
    fig, ax1 = plt.subplots(figsize=figsize)

    ax1.plot(X['z'], X['rms'], color='blue', linewidth=2.5,
             label=r'$\sigma_{x}$')
    ax1.plot(Y['z'], Y['rms'], '--', color='red', linewidth=2.0,
             label=r'$\sigma_{y}$')
    ax1.plot(MAG['z'],
             MAG['profile'] * MAGscale + MAGoffset,
             '-', color='green', linewidth=2.0, label='Magnets')
    ax1.legend()
    ax1.set_ylabel(r'Transverse RMS beam size (mm)', fontsize=22)
    ax1.set_xlabel(r'Distance $z$ (m)', fontsize=22)
    FormatLabelSci()
    plt.tight_layout()
    return fig


# ============================================================================
# Energy Plots
# ============================================================================

def PlotEnergy1plt(X, Y, Z, figsize=(10, 6)):
    """Plot energy, energy spread, and correlated energy spread on dual axes.

    Parameters
    ----------
    X, Y, Z : ndarray
        Emit arrays. Z contains energy info (fields: 'z', 'rmsprime', 'corr', 'avg').
    figsize : tuple
        Figure size.
    """
    fig, ax1 = plt.subplots(figsize=figsize)

    ax1.plot(Z['z'], Z['rmsprime'], color='blue', linewidth=2.0,
             label=r'$\sigma_{E_\mathrm{tot}}$')
    ax1.plot(Z['z'], Z['corr'], color='green', linewidth=2.0,
             label=r'$\sigma_{E_\mathrm{tot}-Cz}$')
    ax1.legend(loc=4)
    ax1.set_ylabel(r'Energy spread $\sigma_E$ (keV)', fontsize=22)
    ax1.set_xlabel(r'Distance $z$ (m)', fontsize=22)
    _add_ymargin_if_flat(ax1, np.concatenate([Z['rmsprime'], Z['corr']]))
    FormatLabelSci()

    ax2 = ax1.twinx()
    ax2.plot(Z['z'], Z['avg'], color='red', linewidth=2.0,
             label=r'$E_\mathrm{kin}$')
    ax2.set_ylabel(r'Kinetic energy (MeV)', fontsize=22, color="red")
    for label in ax2.get_yticklabels():
        label.set_color("red")
    FormatLabelSci()
    plt.tight_layout()
    return fig


# ============================================================================
# Phase Space Density Plots
# ============================================================================

def DensityPlot(X, Y, Nbins, axis=None, cmap=None, figsize=(8, 6)):
    """Create a 2D histogram density plot (log scale).

    Parameters
    ----------
    X, Y : ndarray
        Particle coordinates.
    Nbins : int
        Number of bins per dimension.
    axis : list, optional
        Axis limits [xmin, xmax, ymin, ymax].
    cmap : Colormap, optional
        Default: SLAC-DESY beam map.
    figsize : tuple
        Figure size.

    Returns
    -------
    hist : tuple from hist2d
    """
    if cmap is None:
        cmap = SLAC_DESY_CMAP

    fig = plt.figure(figsize=figsize)
    h = plt.hist2d(X, Y, bins=Nbins, cmap=cmap, norm=LogNorm())

    if axis is not None:
        plt.axis([axis[0], axis[1], axis[2], axis[3]])

    plt.colorbar(label='Counts')
    FormatLabelSci()
    plt.tight_layout()
    return h


def DensityPlot_w_proj(X, Y, Nbins, axis=None, figsize=(8, 6)):
    """Density plot with horizontal projection overlay.

    Parameters
    ----------
    X, Y : ndarray
        Coordinates.
    Nbins : int
        Number of bins.
    axis : list, optional
        Axis limits.
    figsize : tuple
        Figure size.

    Returns
    -------
    hist : tuple from hist2d
    """
    # histogram0 defined at module level

    if axis is not None:
        MyAxis = axis
        ymin, ymax = MyAxis[2], MyAxis[3]
    else:
        ymin, ymax = min(Y), max(Y)

    fig = plt.figure(figsize=figsize)
    h = plt.hist2d(X, Y, bins=Nbins, cmap=SLAC_DESY_CMAP, norm=LogNorm())

    if axis is not None:
        plt.axis(MyAxis)

    plt.colorbar(label='Counts')

    yhist0, xhist0 = histogram0(X, Nbins)
    yhist0 = yhist0 / max(yhist0)
    yhist0 = ymin + (ymax - ymin) * 0.3 * yhist0
    plt.plot(xhist0, yhist0, linewidth=1.5, color='red')

    FormatLabelSci()
    plt.tight_layout()
    return h


def DensityPlot_w_Hproj(X, Y, Nbins, axis=None, figsize=(8, 6)):
    """Density plot with both horizontal and vertical projection overlays.

    Parameters
    ----------
    X, Y : ndarray
        Coordinates.
    Nbins : int
        Number of bins.
    axis : list, optional
        Axis limits.
    figsize : tuple
        Figure size.

    Returns
    -------
    hist : tuple from hist2d
    """
    # histogram0 defined at module level

    if axis is not None:
        MyAxis = axis
        xmin, xmax = MyAxis[0], MyAxis[1]
        ymin, ymax = MyAxis[2], MyAxis[3]
    else:
        xmin, xmax = min(X), max(X)
        ymin, ymax = min(Y), max(Y)

    fig = plt.figure(figsize=figsize)
    h = plt.hist2d(X, Y, bins=Nbins, cmap=SLAC_DESY_CMAP, norm=LogNorm())

    if axis is not None:
        plt.axis(MyAxis)

    plt.colorbar(label='Counts')

    yhist0, xhist0 = histogram0(X, Nbins)
    yhist0 = yhist0 / max(yhist0)
    yhist0 = ymin + (ymax - ymin) * 0.3 * yhist0

    yhist1, xhist1 = histogram0(Y, Nbins)
    yhist1 = yhist1 / max(yhist1)
    yhist1 = xmin + (xmax - xmin) * 0.3 * yhist1

    plt.plot(xhist0, yhist0, linewidth=1.5, color='red')
    plt.plot(yhist1, xhist1, linewidth=1.5, color='red')

    FormatLabelSci()
    plt.tight_layout()
    return h


def DensityplotwProjec2x2(X, Y, Nbins, axis=None, figsize=(12, 10)):
    """2x2 subplot: density + X projection + Y projection.

    Parameters
    ----------
    X, Y : ndarray
        Coordinates.
    Nbins : int
        Number of bins.
    axis : list, optional
        Axis limits [xmin, xmax, ymin, ymax].
    figsize : tuple
        Figure size.
    """
    # histogram0 defined at module level

    fig = plt.figure(figsize=figsize)

    plt.subplot(221)
    if axis is not None:
        MyAxis = axis
        plt.hist2d(X, Y, bins=Nbins, cmap=SLAC_DESY_CMAP, norm=LogNorm())
        plt.axis(MyAxis)
    else:
        MyAxis = [min(X), max(X), min(Y), max(Y)]
        plt.hist2d(X, Y, bins=Nbins, cmap=SLAC_DESY_CMAP, norm=LogNorm())

    plt.colorbar(label='Counts')
    plt.xlabel('X')
    plt.ylabel('Y')
    FormatLabelSci()

    plt.subplot(222)
    yhist0, xhist0 = histogram0(Y, Nbins)
    plt.step(yhist0, xhist0, linewidth=2.5, color='red')
    plt.xlabel('Counts')
    plt.ylabel('Y')
    if axis is not None:
        plt.ylim([MyAxis[2], MyAxis[3]])
    FormatLabelSci()

    plt.subplot(223)
    yhist0, xhist0 = histogram0(X, Nbins)
    plt.step(xhist0, yhist0, linewidth=2.5, color='red')
    plt.xlabel('X')
    plt.ylabel('Counts')
    if axis is not None:
        plt.xlim([MyAxis[0], MyAxis[1]])
    FormatLabelSci()

    plt.tight_layout()
    return fig


# ============================================================================
# Slice Parameter Plots
# ============================================================================

def PlotSliceParameters(sliceMatrix, figsize=(12, 10)):
    """Plot slice analysis results in a 3x2 dashboard.

    Parameters
    ----------
    sliceMatrix : ndarray
        Output from UniformSliceAnalysis (numbins x 14).
        Columns: 0=zcenter, 1=width, 2=Npart, 3=Q, 4=I,
                 5=pav, 6=xcen, 7=ycen, 8=emitnx, 9=emitny,
                 10=dE/E, 11=TB, 12=emitnz, 13=FullB
    figsize : tuple
        Figure size.

    Returns
    -------
    fig : matplotlib Figure
    """
    z = sliceMatrix[:, 0]

    fig, axes = plt.subplots(3, 2, figsize=figsize)

    # Current profile
    axes[0, 0].plot(z, sliceMatrix[:, 4], 'b-', linewidth=2)
    axes[0, 0].set_ylabel('Current (A)')
    axes[0, 0].set_xlabel('z (m)')
    FormatLabelSci()

    # Slice emittance
    axes[0, 1].plot(z, sliceMatrix[:, 8] * 1e6, 'b-', linewidth=2,
                    label=r'$\epsilon_{nx}$')
    axes[0, 1].plot(z, sliceMatrix[:, 9] * 1e6, 'r--', linewidth=2,
                    label=r'$\epsilon_{ny}$')
    axes[0, 1].set_ylabel(r'Slice emittance ($\mu$m)')
    axes[0, 1].set_xlabel('z (m)')
    axes[0, 1].legend()
    FormatLabelSci()

    # Centroid position
    axes[1, 0].plot(z, sliceMatrix[:, 6] * 1e3, 'b-', linewidth=2,
                    label=r'$\langle x \rangle$')
    axes[1, 0].plot(z, sliceMatrix[:, 7] * 1e3, 'r--', linewidth=2,
                    label=r'$\langle y \rangle$')
    axes[1, 0].set_ylabel('Centroid (mm)')
    axes[1, 0].set_xlabel('z (m)')
    axes[1, 0].legend()
    FormatLabelSci()

    # Energy spread
    axes[1, 1].plot(z, sliceMatrix[:, 10] * 100, 'g-', linewidth=2)
    axes[1, 1].set_ylabel(r'$\sigma_E/\langle E \rangle$ (%)')
    axes[1, 1].set_xlabel('z (m)')
    FormatLabelSci()

    # Brightness
    axes[2, 0].plot(z, sliceMatrix[:, 11], 'b-', linewidth=2)
    axes[2, 0].set_ylabel(r'Transverse brightness (A/m$^2$)')
    axes[2, 0].set_xlabel('z (m)')
    FormatLabelSci()

    # Number of particles per slice
    axes[2, 1].plot(z, sliceMatrix[:, 2], 'k-', linewidth=2)
    axes[2, 1].set_ylabel('Particles per slice')
    axes[2, 1].set_xlabel('z (m)')
    FormatLabelSci()

    plt.tight_layout()
    return fig


# ============================================================================
# Bunch Form Factor Plot
# ============================================================================

def PlotBunchFormFactor(k, bff, figsize=(8, 5)):
    """Plot the bunch form factor on log-log axes.

    Parameters
    ----------
    k : ndarray
        Wavenumbers (1/m).
    bff : ndarray
        Form factor values.
    figsize : tuple
        Figure size.
    """
    fig, ax = plt.subplots(figsize=figsize)
    ax.loglog(k, bff, 'b-', linewidth=2)
    ax.set_xlabel(r'Wavenumber $k$ (1/m)')
    ax.set_ylabel(r'Bunch Form Factor $|F(k)|^2$')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    return fig

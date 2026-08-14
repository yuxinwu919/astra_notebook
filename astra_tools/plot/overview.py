"""6D phase-space overview (3x2 grid, postpro style, KDE density).

Column 1: phase spaces x-x', y-y', z-dp/p.
Column 2: spatial correlations x-y, z-x, z-y.
All panels use the unified KDE density engine with 0.5-99.5 percentile
range clipping (outliers can never collapse the display).
"""

from __future__ import annotations

from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from ..analysis.emittance import canonical_divergence
from ..distribution import Distribution
from ._density import clip_percentile, density2d


def plot_overview(
    dist: Distribution,
    bins: int = 120,
    figsize: Tuple[float, float] = (13, 12),
    title: Optional[str] = None,
    bz_on_axis_T: float = 0.0,
    use_weights: bool = False,
) -> Tuple[plt.Figure, np.ndarray]:
    """3x2 overview: phase spaces + spatial correlations.

    Args:
        dist: Distribution.
        bz_on_axis_T: solenoid field at bunch center for canonical x'/y'.
        use_weights: weight densities by macro-particle charge.
    """
    m = dist.active
    p_ref = dist.ref_momentum_eVc
    if p_ref <= 0:
        p_ref = float(np.mean(np.abs(dist.pz[m])))
    w = dist.charge[m] if use_weights else None

    x = dist.x[m] * 1e3
    y = dist.y[m] * 1e3
    z = (dist.z[m] - np.mean(dist.z[m])) * 1e3
    pz = dist.pz[m]
    mean_pz = np.mean(pz)
    dp = (pz - mean_pz) / mean_pz * 100
    ptx = canonical_divergence(dist.px[m], dist.y[m], bz_on_axis_T, +1.0)
    pty = canonical_divergence(dist.py[m], dist.x[m], bz_on_axis_T, -1.0)
    xp = (ptx - np.mean(ptx)) / p_ref * 1e3
    yp = (pty - np.mean(pty)) / p_ref * 1e3

    panels = [
        (0, 0, x, xp, "x-x'", "x [mm]", "x' [mrad]"),
        (0, 1, x, y, "x-y", "x [mm]", "y [mm]"),
        (1, 0, y, yp, "y-y'", "y [mm]", "y' [mrad]"),
        (1, 1, z, x, "z-x", "z [mm]", "x [mm]"),
        (2, 0, z, dp, "z-dp/p", "z [mm] (positive = ahead)", "dp/p [%]"),
        (2, 1, z, y, "z-y", "z [mm]", "y [mm]"),
    ]

    fig, axes = plt.subplots(3, 2, figsize=figsize)
    for row, col, xd, yd, t, xl, yl in panels:
        ax = axes[row, col]
        zz, xe, ye = density2d(xd, yd, bins=bins,
                               range_xy=(clip_percentile(xd), clip_percentile(yd)),
                               weights=w)
        im = ax.pcolormesh(xe, ye, np.where(zz > 0, zz, np.nan).T,
                           cmap="viridis", norm=LogNorm(), shading="auto",
                           rasterized=True)
        fig.colorbar(im, ax=ax)
        ax.axhline(0, color="w", lw=0.6, ls="--", alpha=0.6)
        ax.axvline(0, color="w", lw=0.6, ls="--", alpha=0.6)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_title(t)
    if title:
        fig.suptitle(title, y=1.01)
    fig.tight_layout(pad=2.0, h_pad=1.6, w_pad=1.6)
    return fig, axes


def plot_transverse_profile(
    dist: Distribution,
    bins: int = 140,
    figsize=(6.5, 5),
    title: Optional[str] = None,
) -> plt.Figure:
    """Transverse x-y beam profile (KDE density, log scale)."""
    m = dist.active
    x = dist.x[m] * 1e3
    y = dist.y[m] * 1e3
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    zz, xe, ye = density2d(x, y, bins=bins)
    im = ax.pcolormesh(xe, ye, np.where(zz > 0, zz, np.nan).T,
                       cmap="inferno", norm=LogNorm(), shading="auto",
                       rasterized=True)
    fig.colorbar(im, ax=ax, label="probability density [1/mm^2]")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_aspect("equal")
    ax.set_title(title or "transverse beam profile (x-y)")
    fig.tight_layout()
    return fig

"""Phase space plots: transverse (x-x', y-y') and longitudinal (z-δp/p).

Density mode uses the canonical SLAC-DESY colormap and log-scale normalization
from the embedded ``_plotting`` core.
"""

from __future__ import annotations

from typing import Literal, Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

from ..analysis.emittance import compute_emittance_ellipse_params
from ..distribution import Distribution
from .._plotting.cosmetics import SLAC_DESY_CMAP
from ._artists import draw_emittance_ellipse, draw_reference_lines, add_colorbar
from ._precompute import clip_percentile

Plane = Literal["x", "y", "z"]
Kind = Literal["density", "scatter"]


def plot_phase_space(
    dist: Distribution,
    plane: Plane = "x",
    kind: Kind = "density",
    bins: int = 80,
    ax: Optional[plt.Axes] = None,
    figsize: tuple[float, float] = (6, 5),
    title: Optional[str] = None,
    cmap = None,
    colorbar: bool = True,
    show_ellipse: bool = False,
    use_weights: bool = False,
    **kwargs,
) -> plt.Figure:
    """Plot 2D phase space projection.

    Args:
        dist: Particle distribution.
        plane: Phase plane to plot ('x' → x-x', 'y' → y-y', 'z' → z-δp/p).
        kind: 'density' for 2D histogram, 'scatter' for scatter plot.
        bins: Number of bins for density plot.
        ax: Optional matplotlib Axes. If None, creates a new figure.
        figsize: Figure size (width, height) when creating new figure.
        title: Optional plot title.
        cmap: Colormap name for density plot.
        colorbar: Whether to show a colorbar.
        show_ellipse: Overlay RMS emittance ellipse (only for 'x', 'y' planes).
        use_weights: Weight histogram by macro-particle charge.
        **kwargs: Passed to hist2d or scatter.

    Returns:
        matplotlib Figure.
    """
    mask = dist.active
    pz_abs = np.abs(dist.pz[mask])
    charge_arr = dist.charge[mask] if use_weights else None

    if plane == "x":
        x_data = dist.x[mask] * 1e3  # m → mm
        y_data = (dist.px[mask] - np.mean(dist.px[mask])) / pz_abs * 1e3
        xlabel, ylabel = "x [mm]", "x' [mrad]"
        default_title = "x–x' Phase Space"
    elif plane == "y":
        x_data = dist.y[mask] * 1e3
        y_data = (dist.py[mask] - np.mean(dist.py[mask])) / pz_abs * 1e3
        xlabel, ylabel = "y [mm]", "y' [mrad]"
        default_title = "y–y' Phase Space"
    elif plane == "z":
        x_data = dist.z[mask] * 1e3
        mean_pz = np.mean(dist.pz[mask])
        y_data = (dist.pz[mask] - mean_pz) / mean_pz * 100  # δp/p [%]
        xlabel, ylabel = "z [mm]", "δp/p [%]"
        default_title = "Longitudinal Phase Space"
    else:
        raise ValueError(f"Unknown plane: {plane}. Use 'x', 'y', or 'z'.")

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure

    if kind == "density":
        if cmap is None:
            cmap = SLAC_DESY_CMAP
        vmin_x, vmax_x = clip_percentile(x_data)
        vmin_y, vmax_y = clip_percentile(y_data)
        h = ax.hist2d(
            x_data, y_data, bins=bins, cmap=cmap, norm=LogNorm(),
            range=[[vmin_x, vmax_x], [vmin_y, vmax_y]],
            weights=charge_arr,
            **kwargs,
        )
        if colorbar:
            add_colorbar(fig, ax, h[3])
    else:
        ax.scatter(x_data, y_data, s=1, alpha=0.5, **kwargs)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title or default_title)

    draw_reference_lines(ax)

    if show_ellipse and plane in ("x", "y"):
        x_centered = x_data - np.mean(x_data)
        y_centered = y_data - np.mean(y_data)
        draw_emittance_ellipse(
            ax, x_centered, y_centered,
            n_sigma=1.0, color="red", weights=charge_arr,
        )

    fig.tight_layout()
    return fig


def plot_transverse_phase_space(
    dist: Distribution,
    kind: Kind = "density",
    bins: int = 80,
    figsize: tuple[float, float] = (12, 5),
    title_prefix: str = "",
    cmap = None,
    show_ellipse: bool = False,
    use_weights: bool = False,
) -> plt.Figure:
    """Plot both x-x' and y-y' phase spaces side by side.

    Args:
        dist: Particle distribution.
        kind: 'density' or 'scatter'.
        bins: Number of bins for density plots.
        figsize: Figure size.
        title_prefix: Optional prefix for titles.
        cmap: Colormap name.
        show_ellipse: Overlay emittance ellipses.
        use_weights: Weight by macro-particle charge.

    Returns:
        matplotlib Figure with 2 subplots.
    """
    fig, axes = plt.subplots(1, 2, figsize=figsize)

    for ax, plane in zip(axes, ["x", "y"]):
        plot_phase_space(
            dist, plane=plane, kind=kind, bins=bins,
            ax=ax, cmap=cmap, colorbar=True,
            show_ellipse=show_ellipse,
            use_weights=use_weights,
            title=f"{title_prefix}{'x–x' if plane == 'x' else 'y–y'} Phase Space" if title_prefix else None,
        )

    fig.tight_layout()
    return fig

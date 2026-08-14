"""Six-dimensional phase space overview plot (3×2 grid).

Also usable as a standalone function outside the GUI:
    from beamscope.plot.overview import plot_overview
    fig, axes = plot_overview(dist)
"""

from __future__ import annotations

from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from beamscope.distribution import Distribution
from beamscope.plot._precompute import precompute, clip_percentile, get_overview_panels
from beamscope.plot._artists import add_colorbar, draw_reference_lines


def plot_overview(
    dist: Distribution,
    bins: int = 60,
    cmap: str = "viridis",
    figsize: Tuple[float, float] = (12, 12),
    title: Optional[str] = None,
    use_weights: bool = False,
) -> Tuple[plt.Figure, np.ndarray]:
    """Plot 3×2 six-dimensional phase space overview.

    Column 1 (phase space): x-x', y-y', z-δp/p
    Column 2 (spatial correlations): x-y, z-x, z-y

    All quantities are centroid-subtracted and in display-friendly units
    (mm, mrad, %).

    Args:
        dist: Particle distribution.
        bins: Number of bins for hist2d.
        cmap: Colormap name.
        figsize: Figure size (width, height) in inches.
        title: Optional suptitle.
        use_weights: If True, weight histogram by macro-particle charge.

    Returns:
        (fig, axes) — Figure and 3×2 ndarray of Axes.
    """
    data = precompute(dist)
    panels = get_overview_panels()

    # Charge weights (optional)
    if use_weights:
        mask = dist.active
        weights = dist.charge[mask]
    else:
        weights = None

    fig, axes = plt.subplots(3, 2, figsize=figsize)

    for row, col, x_key, y_key, panel_title, x_label, y_label in panels:
        ax = axes[row, col]
        x_data = data.get(x_key, np.array([]))
        y_data = data.get(y_key, np.array([]))

        if len(x_data) > 0:
            vmin_x, vmax_x = clip_percentile(x_data)
            vmin_y, vmax_y = clip_percentile(y_data)
            h = ax.hist2d(
                x_data, y_data, bins=bins, cmap=cmap,
                range=[[vmin_x, vmax_x], [vmin_y, vmax_y]],
                weights=weights,
            )
            add_colorbar(fig, ax, h[3])
            draw_reference_lines(ax)

        ax.set_xlabel(x_label, fontsize=11)
        ax.set_ylabel(y_label, fontsize=11)
        ax.set_title(panel_title, fontsize=13, fontweight="bold")
        ax.tick_params(labelsize=9)

    if title:
        fig.suptitle(title, fontsize=14, fontweight="bold", y=1.01)

    fig.tight_layout(pad=2.0, h_pad=1.5, w_pad=1.5)
    return fig, axes


def find_panel(row: int, col: int) -> tuple[str, str, str, str, str]:
    """Find the panel definition for a given (row, col).

    Returns (x_key, y_key, title, x_label, y_label).
    """
    for r, c, x_key, y_key, title, x_label, y_label in get_overview_panels():
        if r == row and c == col:
            return x_key, y_key, title, x_label, y_label
    return "x", "xp", "?", "x", "x'"

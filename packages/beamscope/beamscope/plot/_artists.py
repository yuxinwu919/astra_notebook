"""Reusable matplotlib artist elements for beam physics plots.

Provides standardized drawing functions for:
  - RMS emittance ellipse overlay
  - Reference lines (zero-crossing)
  - Statistics text boxes
  - Colorbar helpers
  - Marginal distribution histograms
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Ellipse


# ---------------------------------------------------------------------------
# Emittance ellipse
# ---------------------------------------------------------------------------

def draw_emittance_ellipse(
    ax: plt.Axes,
    x_data: np.ndarray,
    y_data: np.ndarray,
    n_sigma: float = 1.0,
    color: str = "red",
    linewidth: float = 2.0,
    linestyle: str = "-",
    alpha: float = 0.9,
    label: Optional[str] = None,
    annotate: bool = True,
    weights: Optional[np.ndarray] = None,
) -> dict:
    """Draw RMS emittance ellipse on a matplotlib Axes.

    Uses SVD-based computation for numerical stability.

    Args:
        ax: Matplotlib Axes.
        x_data: Centered x-coordinate data (display units).
        y_data: Centered y-coordinate data (display units).
        n_sigma: Number of RMS (1 = RMS ellipse).
        color: Ellipse line color.
        linewidth: Line width.
        linestyle: Line style.
        alpha: Transparency.
        label: Legend label.
        annotate: If True, add text annotation with ε, β, α.
        weights: Optional charge weights for weighted emittance.

    Returns:
        Dict with keys: 'emit', 'beta', 'alpha', 'gamma_t', 'a', 'b', 'theta'.
    """
    from beamscope.analysis.emittance import compute_emittance_ellipse_params

    params = compute_emittance_ellipse_params(
        x_data, y_data, n_sigma=n_sigma, weights=weights,
    )

    if params["a"] <= 0 or params["b"] <= 0:
        return params

    ellipse = Ellipse(
        xy=(0, 0),
        width=2 * params["a"],
        height=2 * params["b"],
        angle=np.degrees(params["theta"]),
        edgecolor=color,
        facecolor="none",
        linewidth=linewidth,
        linestyle=linestyle,
        alpha=alpha,
        label=label or f"{n_sigma}σ RMS",
    )
    ax.add_patch(ellipse)

    if annotate:
        textstr = (
            f"ε={params['emit']:.4f}\n"
            f"β={params['beta']:.2f}\n"
            f"α={params['alpha']:.3f}"
        )
        ax.text(
            0.02, 0.98, textstr, transform=ax.transAxes,
            fontsize=9, verticalalignment="top",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8),
        )

    ax.legend(fontsize=9, loc="upper right")
    return params


# ---------------------------------------------------------------------------
# Reference lines
# ---------------------------------------------------------------------------

def draw_reference_lines(
    ax: plt.Axes,
    x0: float = 0.0,
    y0: float = 0.0,
    color: str = "white",
    linestyle: str = "--",
    linewidth: float = 0.8,
    alpha: float = 0.6,
) -> None:
    """Draw horizontal and vertical reference lines through (x0, y0).

    Args:
        ax: Matplotlib Axes.
        x0: x-intercept.
        y0: y-intercept.
        color: Line color.
        linestyle: Line style.
        linewidth: Line width.
        alpha: Transparency.
    """
    ax.axhline(y0, color=color, linestyle=linestyle, linewidth=linewidth,
               alpha=alpha)
    ax.axvline(x0, color=color, linestyle=linestyle, linewidth=linewidth,
               alpha=alpha)


# ---------------------------------------------------------------------------
# Statistics text box
# ---------------------------------------------------------------------------

def draw_stats_textbox(
    ax: plt.Axes,
    stats: dict,
    position: tuple[float, float] = (0.02, 0.98),
    fontsize: int = 9,
    **kwargs,
) -> None:
    """Draw a statistics summary text box on an Axes.

    Args:
        ax: Matplotlib Axes.
        stats: Dict of label → value pairs.
        position: (x, y) in axes coordinates.
        fontsize: Font size.
        **kwargs: Passed to ax.text().
    """
    lines = [f"{k}: {v}" for k, v in stats.items()]
    textstr = "\n".join(lines)
    ax.text(
        position[0], position[1], textstr,
        transform=ax.transAxes,
        fontsize=fontsize,
        verticalalignment="top",
        bbox=dict(boxstyle="round", facecolor="white", alpha=0.85),
        **kwargs,
    )


# ---------------------------------------------------------------------------
# Colorbar helper
# ---------------------------------------------------------------------------

def add_colorbar(
    fig: plt.Figure,
    ax: plt.Axes,
    mappable,
    label: str = "Count",
    fraction: float = 0.046,
    pad: float = 0.04,
    fontsize: int = 9,
) -> None:
    """Add a standardized colorbar to a plot.

    Args:
        fig: Matplotlib Figure.
        ax: The Axes the mappable belongs to.
        mappable: The image/mesh returned by hist2d, pcolormesh, etc.
        label: Colorbar label.
        fraction: Fraction of axes to use for colorbar.
        pad: Padding between axes and colorbar.
        fontsize: Label font size.
    """
    cbar = fig.colorbar(mappable, ax=ax, fraction=fraction, pad=pad)
    cbar.set_label(label, fontsize=fontsize)
    cbar.ax.tick_params(labelsize=fontsize - 1)


# ---------------------------------------------------------------------------
# Marginal distribution helper
# ---------------------------------------------------------------------------

def draw_marginals(
    fig: plt.Figure,
    x_data: np.ndarray,
    y_data: np.ndarray,
    bins: int = 60,
    color: str = "steelblue",
    alpha: float = 0.7,
    fit_gaussian: bool = True,
) -> tuple[plt.Axes, plt.Axes, plt.Axes, plt.Axes]:
    """Create a jointplot-style layout with marginal histograms.

    Layout:
        ┌──────────┬──────────┐
        │          │  Y marg  │
        │  Main    │          │
        │          ├──────────┤
        │          │ (corner) │
        ├──────────┴──────────┤
        │      X marginal     │
        └─────────────────────┘

    Args:
        fig: Matplotlib Figure.
        x_data: X-axis data.
        y_data: Y-axis data.
        bins: Number of histogram bins.
        color: Fill color.
        alpha: Fill transparency.
        fit_gaussian: Overlay Gaussian fit on marginals.

    Returns:
        (ax_main, ax_top, ax_right, ax_corner) tuple.
    """
    from scipy import stats as sp_stats

    gs = fig.add_gridspec(
        3, 3, hspace=0.05, wspace=0.05,
        width_ratios=[4, 4, 1],
        height_ratios=[1, 4, 4],
    )

    ax_main = fig.add_subplot(gs[1:, :2])
    ax_top = fig.add_subplot(gs[0, :2], sharex=ax_main)
    ax_right = fig.add_subplot(gs[1:, 2], sharey=ax_main)
    ax_corner = fig.add_subplot(gs[0, 2])
    ax_corner.axis("off")

    # X marginal (top)
    ax_top.hist(
        x_data, bins=bins, density=True,
        color=color, alpha=alpha, edgecolor="white",
    )
    mu_x, sig_x = float(np.mean(x_data)), float(np.std(x_data, ddof=1))
    if fit_gaussian:
        xf = np.linspace(x_data.min(), x_data.max(), 200)
        ax_top.plot(xf, sp_stats.norm.pdf(xf, mu_x, sig_x), "r-", lw=1.5)
    ax_top.set_ylabel("Density", fontsize=10)
    ax_top.tick_params(labelbottom=False, labelsize=9)
    ax_top.set_title(f"μ={mu_x:.3f}, σ={sig_x:.3f}", fontsize=10)

    # Y marginal (right)
    ax_right.hist(
        y_data, bins=bins, density=True, orientation="horizontal",
        color=color, alpha=alpha, edgecolor="white",
    )
    mu_y, sig_y = float(np.mean(y_data)), float(np.std(y_data, ddof=1))
    if fit_gaussian:
        yf = np.linspace(y_data.min(), y_data.max(), 200)
        ax_right.plot(sp_stats.norm.pdf(yf, mu_y, sig_y), yf, "r-", lw=1.5)
    ax_right.set_xlabel("Density", fontsize=10)
    ax_right.tick_params(labelleft=False, labelsize=9)
    ax_right.set_title(
        f"μ={mu_y:.3f}, σ={sig_y:.3f}",
        fontsize=10, rotation=-90, x=1.1, y=0.5,
    )

    return ax_main, ax_top, ax_right, ax_corner

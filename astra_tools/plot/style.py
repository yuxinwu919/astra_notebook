"""Matplotlib style presets for modern, publication-quality figures.

No import-time side effects; call set_style() once per notebook.
"""

from __future__ import annotations

from matplotlib import rcParams

# Colorblind-friendly qualitative palette (Paul Tol)
COLORS = ["#0077BB", "#EE7733", "#009988", "#CC3311", "#33BBEE", "#EE3377"]

# Density colormap (default: viridis; alternative SLAC-DESY beam map)
DEFAULT_CMAP = "viridis"


def set_style(
    font_size: int = 12,
    dpi: int = 120,
    fig_width_inches: float = 6.0,
    fig_height_inches: float = 4.5,
) -> None:
    """Apply the astra-notebook figure theme."""
    rcParams.update({
        "figure.dpi": dpi,
        "savefig.dpi": 150,
        "figure.figsize": (fig_width_inches, fig_height_inches),
        "savefig.bbox": "tight",
        "font.size": font_size,
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "axes.labelsize": font_size + 1,
        "axes.titlesize": font_size + 2,
        "axes.linewidth": 1.0,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linestyle": "--",
        "grid.linewidth": 0.6,
        "axes.unicode_minus": False,
        "axes.prop_cycle": __import__("matplotlib").rcsetup.cycler("color", COLORS),
        "xtick.labelsize": font_size - 1,
        "ytick.labelsize": font_size - 1,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "legend.fontsize": font_size - 1,
        "legend.frameon": True,
        "legend.framealpha": 0.9,
        "legend.edgecolor": "0.8",
        "lines.linewidth": 1.6,
        "lines.markersize": 5,
        "image.cmap": DEFAULT_CMAP,
        "text.usetex": False,
    })


def reset_style() -> None:
    """Reset rcParams to matplotlib defaults."""
    import matplotlib.pyplot as plt

    rcParams.update(plt.rcParamsDefault)

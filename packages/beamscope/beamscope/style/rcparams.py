"""Matplotlib rcParams presets for publication-quality figures.

Combines beamscope's comprehensive defaults (DPI, grid, ticks, legend)
with the canonical astra_plotter style.
"""

from __future__ import annotations

import matplotlib.pyplot as plt
from matplotlib import rcParams


def set_publication_style(
    use_tex: bool = False,
    font_size: int = 12,
    dpi: int = 120,
    fig_width_inches: float = 6.0,
    fig_height_inches: float = 4.5,
) -> None:
    """Apply publication-quality matplotlib style.

    Args:
        use_tex: If True, use LaTeX for text rendering (requires LaTeX).
        font_size: Base font size.
        dpi: Figure DPI.
        fig_width_inches: Default figure width.
        fig_height_inches: Default figure height.
    """
    rcParams.update({
        # Figure
        "figure.dpi": dpi,
        "savefig.dpi": 150,
        "figure.figsize": (fig_width_inches, fig_height_inches),
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.05,

        # Font
        "font.size": font_size,
        "font.family": "serif" if use_tex else "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "font.serif": ["Times New Roman", "Computer Modern Roman"],

        # Axes
        "axes.labelsize": font_size + 1,
        "axes.titlesize": font_size + 2,
        "axes.linewidth": 1.0,
        "axes.grid": True,
        "axes.grid.which": "major",
        "grid.alpha": 0.3,
        "grid.linestyle": "--",
        "axes.unicode_minus": False,

        # Ticks
        "xtick.labelsize": font_size - 1,
        "ytick.labelsize": font_size - 1,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.major.size": 5,
        "ytick.major.size": 5,
        "xtick.minor.size": 3,
        "ytick.minor.size": 3,
        "xtick.top": True,
        "ytick.right": True,

        # Legend
        "legend.fontsize": font_size - 1,
        "legend.frameon": True,
        "legend.framealpha": 0.8,
        "legend.edgecolor": "0.5",

        # Lines
        "lines.linewidth": 1.5,
        "lines.markersize": 5,
    })

    if use_tex:
        rcParams.update({
            "text.usetex": True,
            "text.latex.preamble": r"\usepackage{amsmath}\usepackage{siunitx}",
            "font.family": "serif",
        })
    else:
        rcParams.update({"text.usetex": False})

    # Set default colormap from the embedded _plotting core
    try:
        from .._plotting.cosmetics import SLAC_DESY_CMAP
        rcParams["image.cmap"] = SLAC_DESY_CMAP.name
    except ImportError:
        pass


def reset_style() -> None:
    """Reset matplotlib rcParams to defaults."""
    rcParams.update(plt.rcParamsDefault)

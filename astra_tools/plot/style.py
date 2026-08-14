"""Matplotlib style presets for publication-quality figures.

Neutral rcParams suitable for reports and papers; all beam plots should
call set_style() first. Unlike the legacy code, no global side effects
at import time.
"""

from __future__ import annotations

from matplotlib import rcParams


def set_style(
    use_tex: bool = False,
    font_size: int = 12,
    dpi: int = 120,
    fig_width_inches: float = 6.0,
    fig_height_inches: float = 4.5,
) -> None:
    """Apply publication-quality matplotlib rcParams."""
    rcParams.update({
        "figure.dpi": dpi,
        "savefig.dpi": 150,
        "figure.figsize": (fig_width_inches, fig_height_inches),
        "savefig.bbox": "tight",
        "font.size": font_size,
        "font.family": "serif" if use_tex else "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
        "axes.labelsize": font_size + 1,
        "axes.titlesize": font_size + 2,
        "axes.linewidth": 1.0,
        "axes.grid": True,
        "grid.alpha": 0.3,
        "grid.linestyle": "--",
        "axes.unicode_minus": False,
        "xtick.labelsize": font_size - 1,
        "ytick.labelsize": font_size - 1,
        "xtick.direction": "in",
        "ytick.direction": "in",
        "xtick.top": True,
        "ytick.right": True,
        "legend.fontsize": font_size - 1,
        "legend.frameon": True,
        "legend.framealpha": 0.8,
        "lines.linewidth": 1.5,
        "lines.markersize": 5,
    })
    if use_tex:
        rcParams.update({
            "text.usetex": True,
            "text.latex.preamble": r"\usepackage{amsmath}\usepackage{siunitx}",
        })
    else:
        rcParams.update({"text.usetex": False})


def reset_style() -> None:
    """Reset rcParams to matplotlib defaults."""
    import matplotlib.pyplot as plt

    rcParams.update(plt.rcParamsDefault)

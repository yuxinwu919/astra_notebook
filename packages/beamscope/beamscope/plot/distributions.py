"""Projection distribution plots: x, y, z, δp/p histograms with Gaussian fit.

Supports charge-weighted histograms and kernel density estimation.
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats as sp_stats

from ..distribution import Distribution
from ._precompute import precompute


def plot_distributions(
    dist: Distribution,
    bins: int = 80,
    density: bool = True,
    fit_gaussian: bool = True,
    figsize: tuple[float, float] = (14, 4),
    title_prefix: str = "",
    color: str = "steelblue",
    use_weights: bool = False,
    variables: Optional[list[str]] = None,
) -> plt.Figure:
    """Plot projection histograms with optional Gaussian fits.

    Args:
        dist: Particle distribution.
        bins: Number of histogram bins.
        density: If True, normalize to probability density.
        fit_gaussian: If True, overlay Gaussian fit curves.
        figsize: Figure size (width, height).
        title_prefix: Optional prefix for subplot titles.
        color: Histogram fill color.
        use_weights: Weight histogram by macro-particle charge.
        variables: List of variable keys to plot (default: ['x', 'y', 'z']).
            Available: 'x', 'y', 'z', 'xp', 'yp', 'dp', 'E', 'clock', 'r'.

    Returns:
        matplotlib Figure with N subplots.
    """
    data = precompute(dist)

    if variables is None:
        variables = ["x", "y", "z"]

    n = len(variables)
    fig, axes = plt.subplots(1, n, figsize=(5 * n, 4))
    if n == 1:
        axes = [axes]

    if use_weights:
        mask = dist.active
        weights = dist.charge[mask]
    else:
        weights = None

    for ax, key in zip(axes, variables):
        d = data.get(key, np.array([]))
        if len(d) == 0:
            ax.text(0.5, 0.5, "No data", transform=ax.transAxes,
                    ha="center", va="center")
            continue

        ax.hist(
            d, bins=bins, density=density,
            alpha=0.7, color=color, edgecolor="white",
            weights=weights,
        )

        if fit_gaussian:
            if weights is not None:
                mu = float(np.average(d, weights=weights))
                w_sum = np.sum(weights)
                var = np.sum(weights * (d - mu)**2) / (w_sum - np.sum((weights/w_sum)**2))
                sigma = float(np.sqrt(max(var, 0.0)))
            else:
                mu = float(np.mean(d))
                sigma = float(np.std(d, ddof=1))
            x_fit = np.linspace(d.min(), d.max(), 200)
            ax.plot(
                x_fit,
                sp_stats.norm.pdf(x_fit, mu, sigma),
                "r-", linewidth=2,
                label=f"μ={mu:.3f}, σ={sigma:.3f}",
            )
            ax.legend(fontsize=10)

        ax.set_xlabel(f"{key}")
        ax.set_ylabel("Probability Density" if density else "Count")
        ax.set_title(f"{title_prefix}{key} Distribution")

    fig.tight_layout()
    return fig

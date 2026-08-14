"""1D projection distributions with Gaussian fits (postpro style).

x, y, z histograms of the active particles; the z projection is the
longitudinal charge profile. Gaussian overlay uses sample moments.
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from ..distribution import Distribution


def plot_distributions(
    dist: Distribution,
    figsize=(14, 4),
    title_prefix: str = "",
    bins: int = 80,
) -> plt.Figure:
    """x, y, z projection histograms with Gaussian fits.

    All projections are of active particles; x/y/z in mm.
    """
    mask = dist.active
    x = dist.x[mask] * 1e3
    y = dist.y[mask] * 1e3
    z = dist.z[mask] * 1e3

    fig, axes = plt.subplots(1, 3, figsize=figsize)
    for ax, data, label in zip(axes, [x, y, z], ["x", "y", "z"]):
        ax.hist(data, bins=bins, density=True, alpha=0.7,
                color="steelblue", edgecolor="white")
        mu, sigma = float(np.mean(data)), float(np.std(data))
        xs = np.linspace(float(np.min(data)), float(np.max(data)), 200)
        ax.plot(xs, stats.norm.pdf(xs, mu, sigma), "r-", lw=2,
                label="Gauss, sigma=%.4g mm" % sigma)
        ax.set_xlabel(label + " [mm]")
        ax.set_ylabel("probability density")
        ax.set_title((title_prefix + " " + label + " distribution").strip())
        ax.legend(fontsize=8)
    fig.tight_layout()
    return fig


def plot_beam_profile(
    dist: Distribution,
    ax=None,
    figsize=(8, 4),
    bins: int = 100,
    title: Optional[str] = None,
) -> plt.Figure:
    """Transverse x-y beam profile (2D histogram, log scale)."""
    mask = dist.active
    x = dist.x[mask] * 1e3
    y = dist.y[mask] * 1e3
    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure
    h = ax.hist2d(x, y, bins=bins, cmap="inferno")
    fig.colorbar(h[3], ax=ax, label="counts")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_aspect("equal")
    ax.set_title(title or "transverse beam profile (x-y)")
    fig.tight_layout()
    return fig

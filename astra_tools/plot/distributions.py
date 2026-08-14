"""1D projection distributions with Gaussian fits (postpro style).

x, y, z projections of the active particles; the z projection is the
longitudinal charge profile. Gaussian overlays use sample moments;
outliers are folded into the first/last bin (display range clipped to
the 0.2-99.8 percentile).
"""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from ..distribution import Distribution
from ._density import clip_percentile


def _hist_with_outliers(ax, data, bins, color, label_center):
    lo, hi = clip_percentile(data, q=0.2)
    inside = data[(data >= lo) & (data <= hi)]
    n_out = len(data) - len(inside)
    ax.hist(inside, bins=bins, density=True, alpha=0.75,
            color=color, edgecolor="white", linewidth=0.4)
    if n_out:
        ax.text(0.99, 0.98, "%d (%.1f%%) outliers folded" % (n_out, 100 * n_out / len(data)),
                transform=ax.transAxes, ha="right", va="top", fontsize=8, color="0.35")
    return inside, label_center


def plot_distributions(
    dist: Distribution,
    figsize=(14, 4),
    title_prefix: str = "",
    bins: int = 80,
    color: str = "#0077BB",
) -> plt.Figure:
    """x, y, z projection histograms with Gaussian fits (all in mm)."""
    mask = dist.active
    x = dist.x[mask] * 1e3
    y = dist.y[mask] * 1e3
    z = dist.z[mask] * 1e3

    fig, axes = plt.subplots(1, 3, figsize=figsize)
    for ax, data, label in zip(axes, [x, y, z], ["x", "y", "z"]):
        inside, _ = _hist_with_outliers(ax, data, bins, color, label)
        mu, sigma = float(np.mean(inside)), float(np.std(inside))
        xs = np.linspace(float(np.min(inside)), float(np.max(inside)), 200)
        ax.plot(xs, stats.norm.pdf(xs, mu, sigma), color="#CC3311", lw=2,
                label="Gauss: sigma = %.4g mm" % sigma)
        ax.set_xlabel(label + " [mm]")
        ax.set_ylabel("probability density [1/mm]")
        ax.set_title((title_prefix + " " + label + " distribution").strip())
        ax.legend(fontsize=9)
    fig.tight_layout()
    return fig


def plot_energy_distribution(
    dist: Distribution,
    figsize=(6.5, 5),
    bins: int = 80,
    title: str = "energy distribution",
) -> plt.Figure:
    """Kinetic-energy distribution [MeV] with Gaussian fit.

    E_kin = sqrt(pz^2 + m^2 c^4) - m c^2 per particle (relativistic).
    """
    from ..constants import kinetic_energy_from_momentum

    mask = dist.active
    e_kin = kinetic_energy_from_momentum(dist.pz[mask]) * 1e-6  # MeV
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    inside, _ = _hist_with_outliers(ax, e_kin, bins, "#009988", "E")
    mu, sigma = float(np.mean(inside)), float(np.std(inside))
    xs = np.linspace(float(np.min(inside)), float(np.max(inside)), 200)
    ax.plot(xs, stats.norm.pdf(xs, mu, sigma), color="#CC3311", lw=2,
            label="Gauss: sigma = %.4g MeV" % sigma)
    ax.set_xlabel("E_kin [MeV]")
    ax.set_ylabel("probability density [1/MeV]")
    ax.set_title(title)
    ax.legend(fontsize=9)
    fig.tight_layout()
    return fig

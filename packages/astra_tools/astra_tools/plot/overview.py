"""6D phase-space overview (3x2 grid, postpro style).

Column 1: phase spaces x-x', y-y', z-dp/p.
Column 2: spatial correlations x-y, z-x, z-y.
"""

from __future__ import annotations

from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from ..analysis.emittance import canonical_divergence
from ..distribution import Distribution


def _clip(a, q=0.5):
    lo, hi = np.percentile(a, [q, 100.0 - q])
    if hi - lo < 1e-30:
        hi += 1e-30
    return lo, hi


def plot_overview(
    dist: Distribution,
    bins: int = 60,
    figsize: Tuple[float, float] = (12, 12),
    title: Optional[str] = None,
    bz_on_axis_T: float = 0.0,
) -> Tuple[plt.Figure, np.ndarray]:
    """3x2 overview: phase spaces + spatial correlations.

    Args:
        dist: Distribution.
        bz_on_axis_T: solenoid field at bunch center for canonical x'/y'.
    """
    m = dist.active
    x = dist.x[m] * 1e3
    y = dist.y[m] * 1e3
    z = (dist.z[m] - np.mean(dist.z[m])) * 1e3
    pz = dist.pz[m]
    mean_pz = np.mean(pz)
    dp = (pz - mean_pz) / mean_pz * 100
    ptx = canonical_divergence(dist.px[m], dist.y[m], bz_on_axis_T, +1.0)
    pty = canonical_divergence(dist.py[m], dist.x[m], bz_on_axis_T, -1.0)
    xp = (ptx - np.mean(ptx)) / dist.ref_momentum_eVc * 1e3
    yp = (pty - np.mean(pty)) / dist.ref_momentum_eVc * 1e3

    panels = [
        (0, 0, x, xp, "x-x'", "x [mm]", "x' [mrad]"),
        (0, 1, x, y, "x-y", "x [mm]", "y [mm]"),
        (1, 0, y, yp, "y-y'", "y [mm]", "y' [mrad]"),
        (1, 1, z, x, "z-x", "z [mm]", "x [mm]"),
        (2, 0, z, dp, "z-dp/p", "z [mm]", "dp/p [%]"),
        (2, 1, z, y, "z-y", "z [mm]", "y [mm]"),
    ]

    fig, axes = plt.subplots(3, 2, figsize=figsize)
    for row, col, xd, yd, t, xl, yl in panels:
        ax = axes[row, col]
        vx0, vx1 = _clip(xd)
        vy0, vy1 = _clip(yd)
        h = ax.hist2d(xd, yd, bins=bins, cmap="viridis",
                      range=[[vx0, vx1], [vy0, vy1]])
        fig.colorbar(h[3], ax=ax)
        ax.axhline(0, color="w", lw=0.5, ls="--")
        ax.axvline(0, color="w", lw=0.5, ls="--")
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_title(t)
    if title:
        fig.suptitle(title)
    fig.tight_layout(pad=2.0, h_pad=1.5, w_pad=1.5)
    return fig, axes

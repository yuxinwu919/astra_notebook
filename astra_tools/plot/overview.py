"""6D phase-space overview (3x2 grid, postpro style, scatter).

Column 1: phase spaces x-x', y-y', z-dp/p.
Column 2: spatial correlations x-y, z-x, z-y.
所有面板为普通散点图 (确定性子采样), 显示范围 0.5-99.5 百分位
裁剪 (离群点不会压扁主体分布)。
"""

from __future__ import annotations

from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

from ..analysis.emittance import canonical_divergence
from ..distribution import Distribution
from .phase_space import scatter2d


def plot_overview(
    dist: Distribution,
    bins: int = 120,
    figsize: Tuple[float, float] = (13, 12),
    title: Optional[str] = None,
    bz_on_axis_T: float = 0.0,
    use_weights: bool = False,
) -> Tuple[plt.Figure, np.ndarray]:
    """3x2 overview: phase spaces + spatial correlations.

    Args:
        dist: Distribution。
        bz_on_axis_T: solenoid field at bunch center for canonical x'/y'。
        bins/use_weights: 兼容保留 (散点渲染不使用)。
    """
    m = dist.active
    p_ref = dist.ref_momentum_eVc
    if p_ref <= 0:
        p_ref = float(np.mean(np.abs(dist.pz[m])))
    if p_ref <= 0:
        raise ValueError("reference momentum is zero; cannot form x'")

    x = dist.x[m] * 1e3
    y = dist.y[m] * 1e3
    z = (dist.z[m] - np.mean(dist.z[m])) * 1e3
    pz = dist.pz[m]
    mean_pz = np.mean(pz)
    if abs(mean_pz) < 1e-30:
        mean_pz = p_ref
    dp = (pz - mean_pz) / abs(mean_pz) * 100
    ptx = canonical_divergence(dist.px[m], dist.y[m], bz_on_axis_T, +1.0)
    pty = canonical_divergence(dist.py[m], dist.x[m], bz_on_axis_T, -1.0)
    xp = (ptx - np.mean(ptx)) / p_ref * 1e3
    yp = (pty - np.mean(pty)) / p_ref * 1e3

    panels = [
        (0, 0, x, xp, "x-x'", "x [mm]", "x' [mrad]"),
        (0, 1, x, y, "x-y", "x [mm]", "y [mm]"),
        (1, 0, y, yp, "y-y'", "y [mm]", "y' [mrad]"),
        (1, 1, z, x, "z-x", "z [mm]", "x [mm]"),
        (2, 0, z, dp, "z-dp/p", "z [mm] (positive = ahead)", "dp/p [%]"),
        (2, 1, z, y, "z-y", "z [mm]", "y [mm]"),
    ]

    fig, axes = plt.subplots(3, 2, figsize=figsize)
    for row, col, xd, yd, t, xl, yl in panels:
        ax = axes[row, col]
        if len(xd) == 0:
            ax.text(0.5, 0.5, "no active particles",
                    ha="center", va="center", transform=ax.transAxes,
                    color="0.4", fontsize=9)
            ax.set_xlabel(xl)
            ax.set_ylabel(yl)
            ax.set_title(t)
            continue
        scatter2d(ax, xd, yd, max_points=20000)
        ax.axhline(0, color="0.6", lw=0.6, ls="--", alpha=0.6)
        ax.axvline(0, color="0.6", lw=0.6, ls="--", alpha=0.6)
        ax.set_xlabel(xl)
        ax.set_ylabel(yl)
        ax.set_title(t)
    if title:
        fig.suptitle(title, y=1.01)
    fig.tight_layout(pad=2.0, h_pad=1.6, w_pad=1.6)
    return fig, axes


def plot_transverse_profile(
    dist: Distribution,
    bins: int = 140,
    figsize=(6.5, 5),
    title: Optional[str] = None,
) -> plt.Figure:
    """Transverse x-y beam profile (散点图)."""
    m = dist.active
    x = dist.x[m] * 1e3
    y = dist.y[m] * 1e3
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    if len(x) == 0:
        ax.text(0.5, 0.5, "no active particles",
                ha="center", va="center", transform=ax.transAxes,
                color="0.4", fontsize=10)
    else:
        scatter2d(ax, x, y, max_points=30000)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_aspect("equal")
    ax.set_title(title or "transverse beam profile (x-y)")
    fig.tight_layout()
    return fig

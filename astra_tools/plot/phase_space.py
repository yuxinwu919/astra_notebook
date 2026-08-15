"""Phase-space plots (postpro style).

渲染方式: 普通散点图 (点数上限内确定性子采样), 显示范围仍用
0.5-99.5 百分位裁剪, 离群点不会压扁主体分布。

Physics conventions (audited against ASTRA Manual V3.2 and validated
against ASTRA's own Xemit output):
  * x' = p~x / p_ref with canonical momentum p~x = px + c Bz y / 2
    (manual 4.13.1); pass bz_on_axis_T for bunches inside solenoids
  * dp/p from absolute pz (the reader converts the file's relative pz)
  * positive z = ahead of the reference particle (bunch head)
  * normalize=True 除以各自 sigma: 高能束流 (µrad 级发散角) 的
    相空间在原始单位下是一条贴地横线, 归一化后结构清晰可见
"""

from __future__ import annotations

from typing import Optional

import matplotlib.pyplot as plt
import numpy as np

from ..analysis.emittance import canonical_divergence
from ..distribution import Distribution
from ._density import clip_percentile, outside_fraction

PLANES = ("x", "y", "z")


def scatter2d(ax, x, y, max_points=20000, color="#0077BB", s=8.0,
              alpha=0.55, clip_q=0.5):
    """确定性子采样的散点图 + 百分位裁剪范围。

    返回 ((xr, yr), n_total), 供调用方设轴限与离群点注释。
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    if n > max_points:
        idx = np.linspace(0, n - 1, max_points).astype(int)
        x, y = x[idx], y[idx]
    ax.scatter(x, y, s=s, alpha=alpha, color=color,
               edgecolors="none", rasterized=True)
    xr = clip_percentile(x, clip_q)
    yr = clip_percentile(y, clip_q)
    ax.set_xlim(*xr)
    ax.set_ylim(*yr)
    return (xr, yr), n


def plot_phase_space(
    dist: Distribution,
    plane: str = "x",
    ax=None,
    figsize=(6.5, 5),
    title: Optional[str] = None,
    color: str = "#0077BB",
    s: float = 8.0,
    alpha: float = 0.55,
    max_points: int = 20000,
    bz_on_axis_T: float = 0.0,
    clip_q: float = 0.5,
    normalize: bool = False,
    use_weights: bool = False,
) -> plt.Figure:
    """2D 相空间投影散点图。

    Args:
        dist: Distribution。
        plane: 'x' (x-x'), 'y' (y-y') 或 'z' (z-dp/p)。
        max_points: 散点上限 (确定性子采样, 防百万粒子卡顿)。
        bz_on_axis_T: 束心处螺线管轴上场 [T] (正则散角)。
        clip_q: 显示范围百分位裁剪 [%]。
        normalize: 两轴除以各自 sigma (任何能量下结构可见)。
        use_weights: 兼容保留 (散点渲染不使用权重)。
    """
    if plane not in PLANES:
        raise ValueError("plane must be one of " + str(PLANES))

    mask = dist.active

    p_ref = dist.ref_momentum_eVc
    if p_ref <= 0:
        p_ref = float(np.mean(np.abs(dist.pz[mask])))
    if p_ref <= 0:
        raise ValueError("reference momentum is zero; cannot form x'")

    if plane == "x":
        ptx = canonical_divergence(dist.px[mask], dist.y[mask], bz_on_axis_T, +1.0)
        x_data = dist.x[mask] * 1e3
        y_data = (ptx - np.mean(ptx)) / p_ref * 1e3
        xlabel, ylabel = "x [mm]", "x' [mrad]"
        default_title = "x-x' phase space"
    elif plane == "y":
        pty = canonical_divergence(dist.py[mask], dist.x[mask], bz_on_axis_T, -1.0)
        x_data = dist.y[mask] * 1e3
        y_data = (pty - np.mean(pty)) / p_ref * 1e3
        xlabel, ylabel = "y [mm]", "y' [mrad]"
        default_title = "y-y' phase space"
    else:
        x_data = (dist.z[mask] - np.mean(dist.z[mask])) * 1e3
        mean_pz = np.mean(dist.pz[mask])
        if abs(mean_pz) < 1e-30:
            mean_pz = p_ref
        y_data = (dist.pz[mask] - mean_pz) / abs(mean_pz) * 100
        xlabel = "z [mm]  (positive = ahead of reference)"
        ylabel = "dp/p [%]"
        default_title = "longitudinal phase space"

    if normalize:
        sx = float(np.std(x_data))
        sy = float(np.std(y_data))
        if sx > 0:
            x_data = (x_data - np.mean(x_data)) / sx
        if sy > 0:
            y_data = (y_data - np.mean(y_data)) / sy
        if plane in ("x", "y"):
            xlabel, ylabel = plane + "/σ" + plane, plane + "'/σ" + plane + "'"
        else:
            xlabel, ylabel = "z/σz", "(dp/p)/σ"
        default_title += " (normalized)"

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure

    (xr, yr), n_total = scatter2d(ax, x_data, y_data, max_points=max_points,
                                  color=color, s=s, alpha=alpha,
                                  clip_q=clip_q)
    ax.axhline(0, color="0.6", lw=0.6, ls="--", alpha=0.7)
    ax.axvline(0, color="0.6", lw=0.6, ls="--", alpha=0.7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title or default_title)

    f_out = outside_fraction(x_data, *xr) + outside_fraction(y_data, *yr)
    notes = []
    if n_total > max_points:
        notes.append("showing %d / %d points" % (max_points, n_total))
    if f_out > 0.001:
        notes.append("range clip: %.1f%% outside" % (100 * f_out))
    if notes:
        ax.text(0.99, 0.02, "; ".join(notes), transform=ax.transAxes,
                ha="right", va="bottom", fontsize=8, color="0.35")

    fig.tight_layout()
    return fig


def plot_transverse_phase_space(
    dist: Distribution,
    figsize=(13, 5),
    title_prefix: str = "",
    **kwargs,
) -> plt.Figure:
    """Both transverse phase spaces side by side (x-x', y-y')."""
    fig, axes = plt.subplots(1, 2, figsize=figsize)
    plot_phase_space(dist, plane="x", ax=axes[0],
                     title=(title_prefix + " x-x'").strip(), **kwargs)
    plot_phase_space(dist, plane="y", ax=axes[1],
                     title=(title_prefix + " y-y'").strip(), **kwargs)
    return fig

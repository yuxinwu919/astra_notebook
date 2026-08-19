"""任意参数对相空间 (postpro 5.6.2 菜单 2).

从分布中任选两个参数作横/纵坐标散点; 支持加投影 (add_proj)、
减线性相关 (subtract_corr)、状态着色 (color_by_status)。

可用参数 (active 粒子):
  x, y, z  [mm], px/py/pz [MeV/c], clock [ps], t [ps],
  x', y' [mrad] (正则散角), dp/p [%], E_kin [MeV]
"""

from __future__ import annotations

from typing import Optional

import numpy as np
import matplotlib.pyplot as plt

from ..analysis.emittance import canonical_divergence, canonical_signs
from ..analysis.time import bunch_time
from ..constants import kinetic_energy_from_momentum_vector
from ..distribution import Distribution
from .phase_space import scatter2d


def param_columns(dist: Distribution, bz_on_axis_T: float = 0.0,
                  mask=None) -> dict:
    """可用参数列: {name: (values, unit, axis_label)}.

    values 为显示单位 (mm / mrad / % / MeV / ps)。
    mask (R2-2-2): 抽列用的粒子掩码, 默认 active。状态着色时由
    plot_arbitrary 传入 status>=-6 的同一 mask — mask 与数据列
    必须同源, 否则 (passive/lost 粒子存在时) 布尔索引长度不匹配。
    """
    if mask is None:
        mask = dist.active
    p_ref = dist.ref_momentum_or_mean()
    s_can = canonical_signs(dist)[mask]   # 种类感知符号 (2026-08 F4)
    ptx = canonical_divergence(dist.px[mask], dist.y[mask], bz_on_axis_T, s_can)
    pty = canonical_divergence(dist.py[mask], dist.x[mask], bz_on_axis_T, -s_can)
    e = kinetic_energy_from_momentum_vector(
        dist.px[mask], dist.py[mask], dist.pz[mask]) * 1e-6  # MeV (全动量, 2026-08 P2-2)
    dp = (dist.pz[mask] - np.mean(dist.pz[mask])) / abs(np.mean(dist.pz[mask])) * 100
    return {
        "x":     (dist.x[mask] * 1e3, "mm", "x [mm]"),
        "y":     (dist.y[mask] * 1e3, "mm", "y [mm]"),
        "z":     (dist.z[mask] * 1e3, "mm", "z [mm]"),
        "px":    (dist.px[mask] * 1e-6, "MeV/c", "px [MeV/c]"),
        "py":    (dist.py[mask] * 1e-6, "MeV/c", "py [MeV/c]"),
        "pz":    (dist.pz[mask] * 1e-6, "MeV/c", "pz [MeV/c]"),
        "clock": (dist.clock[mask] * 1e12, "ps", "clock [ps]"),
        "t":     (bunch_time(dist)[mask] * 1e12, "ps", "t [ps]"),
        "xp":    ((ptx - np.mean(ptx)) / p_ref * 1e3, "mrad", "x' [mrad]"),
        "yp":    ((pty - np.mean(pty)) / p_ref * 1e3, "mrad", "y' [mrad]"),
        "dp/p":  (dp, "%", "dp/p [%]"),
        "E_kin": (e, "MeV", "E_kin [MeV]"),
    }


def subtract_linear_corr(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """减去 y 对 x 的线性相关 (最小二乘), 返回残差."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if len(x) < 3:
        return y - np.mean(y)
    b, a = np.polyfit(x, y, 1)
    return y - (a + b * x)


def plot_arbitrary(
    dist: Distribution,
    x_param: str,
    y_param: str,
    subtract_corr: bool = False,
    add_proj: bool = False,
    color_by_status: bool = False,
    bz_on_axis_T: float = 0.0,
    figsize=(7, 6),
    title: Optional[str] = None,
    clip_q: float = 0.5,
) -> plt.Figure:
    """任意两参数相空间散点图 (postpro 5.6.2).

    Args:
        dist: 分布。
        x_param / y_param: param_columns() 的键。
        subtract_corr: 减线性相关 (拟合 v=a+b*u 画残差)。
        add_proj: 叠加两轴边缘直方图投影。
        color_by_status: 按手册 Table 6 状态着色。
        bz_on_axis_T: 螺线管轴上场 (正则散角)。
    """
    # R2-2-2: mask 与数据列同源 — 先定 mask, 再让 param_columns
    # 按同一 mask 抽列 (旧实现列来自 active, mask 作用全粒子,
    # passive/lost 存在时布尔索引长度不匹配 ValueError)。
    if color_by_status:
        mask = dist.status >= -6
    else:
        mask = dist.active
    cols = param_columns(dist, bz_on_axis_T, mask=mask)
    if x_param not in cols or y_param not in cols:
        raise KeyError("未知参数 %r / %r, 可用: %s"
                       % (x_param, y_param, ", ".join(cols)))
    x, _, xlab = cols[x_param]
    y, _, ylab = cols[y_param]
    if subtract_corr:
        y = subtract_linear_corr(x, y)
        ylab += " (lin. corr. removed)"

    fig, ax = plt.subplots(figsize=figsize)
    scatter2d(ax, x, y, clip_q=clip_q,
              status=dist.status[mask] if color_by_status else None)
    ax.axhline(0, color="0.6", lw=0.6, ls="--", alpha=0.7)
    ax.axvline(0, color="0.6", lw=0.6, ls="--", alpha=0.7)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title or "%s vs %s" % (y_param, x_param))
    if color_by_status:
        ax.legend(fontsize=8, loc="best", markerscale=1.5)
    if add_proj:
        _add_projections(fig, ax, x, y)
        fig.subplots_adjust(left=0.13, right=0.86, top=0.86, bottom=0.12)
    else:
        fig.tight_layout()
    return fig


def _add_projections(fig, ax, x, y):
    """顶部 + 右侧边缘直方图 (postpro 'add projections')."""
    def _inset(pos):
        return fig.add_axes(pos, frameon=False)

    # 顶部投影 (x)
    axt = _inset([0.15, 0.86, 0.7, 0.09])
    axt.hist(x, bins=40, color="#0077BB", alpha=0.6)
    axt.set_xticks([]); axt.set_yticks([])
    # 右侧投影 (y)
    axr = _inset([0.88, 0.14, 0.09, 0.7])
    axr.hist(y, bins=40, orientation="horizontal", color="#0077BB", alpha=0.6)
    axr.set_xticks([]); axr.set_yticks([])


class OverlayManager:
    """相空间叠加管理器 (postpro 5.6.2 'save & read overlay').

    保存若干分布的散点数据, 叠加画到同一轴; 叠加层用
    反向色 (默认) 或黑色 (Plot_mode=1 风格), 以与首层区分。
    """

    _PALETTE = ["#0077BB", "#CC3311", "#009988", "#EE7733", "#882255"]

    def __init__(self):
        self._layers = []   # [(x, y, xlab, ylab, label)]

    def add(self, dist: Distribution, x_param: str, y_param: str,
            bz_on_axis_T: float = 0.0, label: Optional[str] = None) -> int:
        """保存当前分布的一层散点数据; 返回层序号."""
        cols = param_columns(dist, bz_on_axis_T)
        x, _, xlab = cols[x_param]
        y, _, ylab = cols[y_param]
        self._layers.append((np.asarray(x, float), np.asarray(y, float),
                             xlab, ylab, label or ("layer %d" % len(self._layers))))
        return len(self._layers) - 1

    @property
    def count(self) -> int:
        return len(self._layers)

    def clear(self) -> None:
        self._layers.clear()

    def plot(self, figsize=(7, 6), title=None, clip_q=0.5):
        """把全部层画到同一轴 (首层实线色, 叠加层按反色/黑)."""
        if not self._layers:
            raise ValueError("overlay manager is empty")
        fig, ax = plt.subplots(figsize=figsize)
        allx, ally = [], []
        for i, (x, y, xlab, ylab, lab) in enumerate(self._layers):
            if i == 0:
                c, mk = self._PALETTE[0], "o"
            else:
                # 叠加层: 与首层颜色相反 + 空心标记 (区分度)
                c, mk = (self._PALETTE[i % len(self._PALETTE)] if i > 1
                         else "#BBBBBB"), ["o", "^", "s", "D", "v"][i % 5]
            ax.scatter(x, y, s=8, alpha=0.5, color=c, marker=mk,
                       edgecolors="none", rasterized=True, label=lab)
            allx.append(x); ally.append(y)
        xr = np.percentile(np.concatenate(allx), [0.5, 99.5])
        yr = np.percentile(np.concatenate(ally), [0.5, 99.5])
        ax.set_xlim(*xr); ax.set_ylim(*yr)
        ax.set_xlabel(self._layers[0][2])
        ax.set_ylabel(self._layers[0][3])
        ax.set_title(title or "overlaid phase spaces (n=%d)" % len(self._layers))
        ax.legend(fontsize=8)
        fig.tight_layout()
        return fig

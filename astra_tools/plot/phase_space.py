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
from ..analysis.time import bunch_time
from ..distribution import Distribution
from ._density import clip_percentile, outside_fraction

PLANES = ("x", "y", "z", "t")

# 手册 Table 6: postpro 绘图的状态颜色/符号编码.
# (mask_pred, label, color, marker)
STATUS_SPEC = [
    (">5", "secondary", "#2ca02c", "o"),
    ("2,3,5", "normal", "#111111", "*"),
    ("4", "marked", "#d62728", "*"),
    ("0,1", "passive", "#2ca02c", "+"),
    ("-1..-6", "cathode", "#8c564b", "x"),
    ("-12..-25", "lost aperture", "#ff0000", "."),
    ("-26..-30", "lost", "#1f77b4", "o"),
    ("<=-30", "lost deep", "#1f77b4", "*"),
]


def status_mask(s: np.ndarray, pred: str) -> np.ndarray:
    """按 Table 6 谓词取状态掩码."""
    s = np.asarray(s)
    if pred == ">5":
        return s > 5
    if pred == "2,3,5":
        return (s == 2) | (s == 3) | (s == 5)
    if pred == "4":
        return s == 4
    if pred == "0,1":
        return (s == 0) | (s == 1)
    if pred == "-1..-6":
        return (s >= -6) & (s <= -1)
    if pred == "-12..-25":
        return (s >= -25) & (s <= -12)
    if pred == "-26..-30":
        return (s >= -30) & (s <= -26)
    if pred == "<=-30":
        return s <= -30
    raise ValueError("unknown status predicate: %r" % pred)


def scatter2d(ax, x, y, max_points=20000, color="#0077BB", s=8.0,
              alpha=0.55, clip_q=0.5, status=None):
    """确定性子采样的散点图 + 百分位裁剪范围.

    返回 ((xr, yr), n_total), 供调用方设轴限与离群点注释。
    status 非 None 时按手册 Table 6 分组着色 (每组合适的颜色/marker,
    带 label, 供 ax.legend()); 未匹配状态的粒子用给定 color。
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = len(x)
    idx = np.arange(n)
    color_arr = isinstance(color, (list, tuple, np.ndarray))
    if n > max_points:
        idx = np.linspace(0, n - 1, max_points).astype(int)
        x, y = x[idx], y[idx]
        if color_arr:
            color = np.asarray(color)[idx]
    if status is not None:
        st = np.asarray(status)
        if len(st) == n:
            st = st[idx]
        matched = np.zeros(n, dtype=bool)
        for pred, _label, c, mk in STATUS_SPEC:
            mask = status_mask(st, pred)
            if np.any(mask):
                kw = dict(s=s, alpha=alpha, color=c, marker=mk,
                          rasterized=True)
                if mk not in ("+", "x"):   # unfilled marker 无边缘色
                    kw["edgecolors"] = "none"
                ax.scatter(x[mask], y[mask], label=_label, **kw)
            matched |= mask
        if not np.all(matched):
            ax.scatter(x[~matched], y[~matched], s=s, alpha=alpha,
                       color=color, edgecolors="none", rasterized=True)
    else:
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
    color_by_status: bool = False,
    colors=None,
) -> plt.Figure:
    """2D 相空间投影散点图。

    Args:
        dist: Distribution。
        plane: 'x' (x-x'), 'y' (y-y'), 'z' (z-dp/p) 或 't' (t-dp/p,
            postpro 5.6.1 项 3 时间坐标版)。
        max_points: 散点上限 (确定性子采样, 防百万粒子卡顿)。
        bz_on_axis_T: 束心处螺线管轴上场 [T] (正则散角)。
        clip_q: 显示范围百分位裁剪 [%]。
        normalize: 两轴除以各自 sigma (任何能量下结构可见)。
        use_weights: 兼容保留 (散点渲染不使用权重)。
        color_by_status: 按手册 Table 6 状态分组着色 (含 passive/lost)。
        colors: 可选逐粒子颜色数组 (长度 = 分布粒子数, 如 Plot_steering.par
            CP_ind 的 cp_index_colors); 提供时覆盖单色 (优先于 color_by_status
            之外的默认色)。
    """
    if plane not in PLANES:
        raise ValueError("plane must be one of " + str(PLANES))

    # 手册 5.6: 默认绘图条件 status >= -6 (含 passive/cathode);
    # 常规图仍只用 active (status>1, 统计口径), 状态着色图按手册标准。
    mask = dist.status >= -6 if color_by_status else dist.active

    if not np.any(mask):
        raise ValueError(
            "no active particles (status>1): cannot plot phase space "
            "(batch-2 fix: empty-bunch fallback used to produce nan, "
            "and nan<=0 is False so no error was raised)")
    p_ref = dist.ref_momentum_or_mean()

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
    elif plane == "t":
        # postpro 5.6.1 项 3: 时间坐标版纵向相空间
        x_data = bunch_time(dist)[mask] * 1e12   # ps
        mean_pz = np.mean(dist.pz[mask])
        if abs(mean_pz) < 1e-30:
            mean_pz = p_ref
        y_data = (dist.pz[mask] - mean_pz) / abs(mean_pz) * 100
        xlabel = "t [ps]  (0 = bunch centre)"
        ylabel = "dp/p [%]"
        default_title = "longitudinal phase space (time)"
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
        elif plane == "t":
            xlabel, ylabel = "t/σt", "(dp/p)/σ"
        else:
            xlabel, ylabel = "z/σz", "(dp/p)/σ"
        default_title += " (normalized)"

    if ax is None:
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    else:
        fig = ax.figure

    (xr, yr), n_total = scatter2d(ax, x_data, y_data, max_points=max_points,
                                  color=colors[mask] if colors is not None
                                  else color,
                                  s=s, alpha=alpha,
                                  clip_q=clip_q,
                                  status=dist.status[mask]
                                  if color_by_status else None)
    ax.axhline(0, color="0.6", lw=0.6, ls="--", alpha=0.7)
    ax.axvline(0, color="0.6", lw=0.6, ls="--", alpha=0.7)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title or default_title)
    if color_by_status:
        ax.legend(fontsize=8, loc="best", markerscale=1.5)

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

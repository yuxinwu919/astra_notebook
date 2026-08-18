"""Field-map plots (fieldplot replacement).

Cavity: on-axis Ez(z) and 2D (r, z) maps of Ez, Er, Bphi reconstructed
from the axis field via the off-axis expansion (manual chapter 8).
Solenoid: on-axis Bz(z) and 2D (r, z) map of Bz with Br arrows.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import Optional, Union

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (注册 3d 投影)

from ..io.field_map import (CavityField, SolenoidField, FieldMap3D,
                            read_3d_field_map_components)


def plot_cavity_field(
    field: CavityField,
    omega: float = 0.0,
    rmax: Optional[float] = None,
    figsize=(14, 8),
    title: Optional[str] = None,
    maxE_MVpm: Optional[float] = None,
) -> plt.Figure:
    """Cavity field: Ez(r,z), Er(r,z), Bphi(r,z) maps + on-axis Ez(z).

    Args:
        field: CavityField (on-axis table).
        omega: RF angular frequency [rad/s] for Bphi; 0 = static.
        rmax: max radius for the map [m]; default 10% of cavity length.
        maxE_MVpm: scale the field to this peak [MV/m] (like ASTRA's
            MaxE(n)); None = raw file values.
    """
    # 批 3: ez0 为原始任意单位 (手册 6.9); 给 maxE 时按峰值缩放并
    # 以 MV/m 显示, 否则按任意单位显示 (不再假装 V/m)。
    peak = float(np.max(np.abs(field.ez0)))
    if maxE_MVpm is not None:
        scale = maxE_MVpm / peak if peak > 0 else 0.0
        ez_unit, b_unit = "MV/m", "T"
    else:
        scale = 1.0
        ez_unit, b_unit = "arb. units", "arb. units"

    z0, z1 = float(field.z.min()), float(field.z.max())
    if rmax is None:
        rmax = 0.1 * (z1 - z0)
    nz, nr = 400, 60
    zz = np.linspace(z0, z1, nz)
    rr = np.linspace(0, rmax, nr)
    ZZ, RR = np.meshgrid(zz, rr)
    ez, er, bphi = field.field_at(RR, ZZ, omega)
    ez = ez * scale
    er = er * scale
    bphi = bphi * scale

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    maps = [(ez, "Ez [%s]" % ez_unit), (er, "Er [%s]" % ez_unit),
            (bphi, "Bphi [%s]" % b_unit)]
    for ax, (m, label) in zip(axes.flat[:3], maps):
        vmax = float(np.max(np.abs(m)))
        if vmax == 0:
            vmax = 1.0
        im = ax.pcolormesh(ZZ * 1e3, RR * 1e3, m, cmap="RdBu_r",
                           vmin=-vmax, vmax=vmax, shading="auto")
        fig.colorbar(im, ax=ax, label=label)
        ax.set_xlabel("z [mm]")
        ax.set_ylabel("r [mm]")

    ax = axes[1, 1]
    ax.plot(field.z * 1e3, field.ez0 * scale, color="C0")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("Ez on axis [%s]" % ez_unit)
    ax.set_title("on-axis field")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_solenoid_field(
    field: SolenoidField,
    rmax: Optional[float] = None,
    figsize=(12, 4),
    title: Optional[str] = None,
) -> plt.Figure:
    """Solenoid: Bz(r,z) map and on-axis Bz(z).

    Args:
        field: SolenoidField, already scaled to Tesla (use .scaled(maxB)).
    """
    z0, z1 = float(field.z.min()), float(field.z.max())
    if rmax is None:
        rmax = 0.05 * (z1 - z0)
    nz, nr = 400, 40
    zz = np.linspace(z0, z1, nz)
    rr = np.linspace(0, rmax, nr)
    ZZ, RR = np.meshgrid(zz, rr)
    br, bz = field.field_at(RR, ZZ)

    fig, axes = plt.subplots(1, 2, figsize=figsize)
    ax = axes[0]
    vmax = float(np.max(np.abs(bz))) or 1.0
    im = ax.pcolormesh(ZZ * 1e3, RR * 1e3, bz, cmap="RdBu_r",
                       vmin=-vmax, vmax=vmax, shading="auto")
    fig.colorbar(im, ax=ax, label="Bz [T]")
    # Br arrows on a coarse grid
    stride = 12
    ax.quiver(ZZ[::stride, ::stride] * 1e3, RR[::stride, ::stride] * 1e3,
              np.zeros_like(ZZ[::stride, ::stride]), br[::stride, ::stride],
              color="k", scale=abs(br).max() * 2 or 1.0, width=0.002)
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("r [mm]")
    ax.set_title("Bz(r,z)")

    ax = axes[1]
    ax.plot(field.z * 1e3, field.bz0, color="C0")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("Bz on axis [T]")
    ax.set_title("on-axis field")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_solenoid_components(
    field,
    r_probe: Optional[float] = None,
    figsize=(12, 8),
    title: Optional[str] = None,
) -> plt.Figure:
    """螺线管 next page (fieldplot 菜单 4): 单独 Bz/Br + 场展开半径 R3rd.

    2x2: [轴上 Bz(z), 单独 Br(r_probe, z) 径向梯度],
         [Bz(r,z) 二维图, R3rd(z) 静磁半径 (手册 8 章)]。
    r_probe: 径向探针位置 [m] (默认 5% 纵向跨度)。
    """
    from ..io.field_map import SolenoidField
    if not isinstance(field, SolenoidField):
        raise TypeError("plot_solenoid_components 需要 SolenoidField")
    if r_probe is None:
        span = float(field.z.max() - field.z.min())
        r_probe = 0.05 * span

    z_ax = np.linspace(field.z.min(), field.z.max(), 400)
    br_probe, bz_ax = field.field_at(np.full_like(z_ax, r_probe), z_ax)

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    # 单独 Bz (轴上)
    axes[0, 0].plot(field.z * 1e3, field.bz0, color="C0")
    axes[0, 0].set_xlabel("z [mm]")
    axes[0, 0].set_ylabel("Bz on axis [T]")
    axes[0, 0].set_title("Bz alone")
    # 单独 Br (径向探针)
    axes[0, 1].plot(z_ax * 1e3, br_probe, color="C1")
    axes[0, 1].set_xlabel("z [mm]")
    axes[0, 1].set_ylabel("Br(r=%.1f mm) [T]" % (r_probe * 1e3))
    axes[0, 1].set_title("Br alone (radial gradient)")
    # Bz(r,z) 二维
    rr = np.linspace(0, 3 * r_probe, 40)
    ZZ, RR = np.meshgrid(z_ax, rr)
    br2, bz2 = field.field_at(RR, ZZ)
    vmax = float(np.max(np.abs(bz2)))
    if vmax == 0:
        vmax = 1.0
    im = axes[1, 0].pcolormesh(ZZ * 1e3, RR * 1e3, bz2, cmap="RdBu_r",
                               vmin=-vmax, vmax=vmax, shading="auto")
    fig.colorbar(im, ax=axes[1, 0], label="Bz [T]")
    axes[1, 0].set_xlabel("z [mm]")
    axes[1, 0].set_ylabel("r [mm]")
    axes[1, 0].set_title("Bz(r,z)")
    # R3rd (静磁, 轻度平滑抑制 Bz''' 噪声)
    r3 = np.asarray(field.expansion_radius(smooth_window=5), float)
    ok = np.isfinite(r3) & (r3 >= 0)
    axes[1, 1].plot(field.z[ok] * 1e3, r3[ok] * 1e3, color="C2")
    axes[1, 1].set_xlabel("z [mm]")
    axes[1, 1].set_ylabel(r"$R_{3rd}$ [mm]")
    # R3rd 在 Bz''' 过零处发散: 裁剪 Y 上限到 99 分位 (手册 8 章噪声敏感)
    if ok.any():
        axes[1, 1].set_ylim(0, float(np.percentile(r3[ok], 99)) * 1e3)
    axes[1, 1].set_title("field expansion radius (static)")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_te_field(
    field,
    omega: float = 0.0,
    rmax: Optional[float] = None,
    figsize=(14, 8),
    title: Optional[str] = None,
    max_field: Optional[float] = None,
) -> plt.Figure:
    """TE 模场 (fieldplot 菜单 1 项 2, 手册 8 章 TE 展开).

    Bz(r,z) / Br(r,z) / Eφ(r,z) 二维图 + 轴上 Bz。
    field: TEField (手册 6.9: 文件名 'TE_' 前缀, 表存轴上纵向磁场)。
    omega: RF 角频率 [rad/s] (Eφ 需 omega; 0 = 静态, Eφ=0)。
    max_field: 按峰值缩放 (TE 模的 MaxE 指轴上 Bz); None = 任意单位。
    """
    from ..io.field_map import TEField
    if not isinstance(field, TEField):
        raise TypeError("plot_te_field 需要 TEField")
    peak = float(np.max(np.abs(field.bz0)))
    if max_field is not None:
        scale = max_field / peak if peak > 0 else 0.0
        b_unit, e_unit = "T", "V/m"
    else:
        scale = 1.0
        b_unit, e_unit = "arb. units", "arb. units"

    z0, z1 = float(field.z.min()), float(field.z.max())
    if rmax is None:
        rmax = 0.1 * (z1 - z0)
    zz = np.linspace(z0, z1, 400)
    rr = np.linspace(0, rmax, 60)
    ZZ, RR = np.meshgrid(zz, rr)
    bz, br, ephi = field.field_at(RR, ZZ, omega)
    bz = bz * scale
    br = br * scale
    ephi = ephi * scale

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    maps = [(bz, "Bz [%s]" % b_unit), (br, "Br [%s]" % b_unit),
            (ephi, "Ephi [%s]" % e_unit)]
    for ax, (m, label) in zip(axes.flat[:3], maps):
        vmax = float(np.max(np.abs(m)))
        if vmax == 0:
            vmax = 1.0
        im = ax.pcolormesh(ZZ * 1e3, RR * 1e3, m, cmap="RdBu_r",
                           vmin=-vmax, vmax=vmax, shading="auto")
        fig.colorbar(im, ax=ax, label=label)
        ax.set_xlabel("z [mm]")
        ax.set_ylabel("r [mm]")

    ax = axes[1, 1]
    ax.plot(field.z * 1e3, field.bz0 * scale, color="C0")
    ax.set_xlabel("z [mm]")
    ax.set_ylabel("Bz on axis [%s]" % b_unit)
    ax.set_title("on-axis field")
    if title:
        fig.suptitle(title)
    fig.tight_layout()
    return fig


def plot_field_expansion_radius(
    field,
    omega: float = 0.0,
    figsize=(8, 4),
    title: Optional[str] = None,
) -> plt.Figure:
    """场展开半径 R3rd vs z (手册 8 章诊断量, 对场表数值噪声敏感).

    支持 CavityField (TM: R3rd_Er / R3rd_Bφ)、TEField (R3rd_Eφ / R3rd_Br)、
    SolenoidField (静磁 R3rd = sqrt(|0.08 Bz'|/|Bz'''|))。
    """
    from ..io.field_map import SolenoidField, TEField
    z = np.asarray(field.z, dtype=float)
    if isinstance(field, SolenoidField):
        specs = [(np.asarray(field.expansion_radius(smooth_window=5), float),
                  r"$R_{3rd}$ (solenoid)")]
    elif isinstance(field, TEField):
        r_ep, r_br = field.expansion_radius(omega, smooth_window=5)
        specs = [(np.asarray(r_ep, float), r"$R_{3rd}^{E_\phi}$"),
                 (np.asarray(r_br, float), r"$R_{3rd}^{B_r}$")]
    else:  # CavityField (TM)
        r_er, r_bp = field.expansion_radius(omega, smooth_window=5)
        specs = [(np.asarray(r_er, float), r"$R_{3rd}^{E_r}$"),
                 (np.asarray(r_bp, float), r"$R_{3rd}^{B_\phi}$")]
    fig, ax = plt.subplots(figsize=figsize)
    for r, lab in specs:
        r = np.asarray(r, float)
        r = np.clip(r, 0, None)          # 数值噪声/负值压到 >=0
        ok = np.isfinite(r)              # 0/0 或分母过零处标 NaN, 跳过
        if ok.any():
            ax.plot(z[ok] * 1e3, r[ok] * 1e3, label=lab)
    ax.set_xlabel("z [mm]")
    ax.set_ylabel(r"$R_{3rd}$ [mm]")
    # R3rd 在导数过零处发散 (手册 8 章: 对数值噪声敏感):
    # 裁剪 Y 上限到各曲线 99.5 分位, 避免尖峰压扁物理区间
    allvals = np.concatenate([np.asarray(r_, float)[np.isfinite(r_)
                                                  & (np.asarray(r_) >= 0)]
                              for r_, _ in specs])
    if allvals.size:
        ax.set_ylim(0, float(np.percentile(allvals, 99.5)) * 1e3)
    ax.set_title(title or "field expansion radius R3rd (manual ch. 8)")
    ax.legend(fontsize=9)
    fig.tight_layout()
    return fig


# ---------------------------------------------------------------------------
# 3D 场图: 矢量剖面 / 2D 等值线 / 3D 等值线 / 统一分发器
# ---------------------------------------------------------------------------

_AXIS_IDX = {"x": 0, "y": 1, "z": 2}

# 显示平面 -> (固定轴下标, 横轴, 纵轴, 横轴分量, 纵轴分量, x标签, y标签)
# 约定: xz/yz 平面中 z 在横轴 (加速器物理习惯)。
_PLANE_SPEC = {
    "xy": (2, "x", "y", "fx", "fy", "x [mm]", "y [mm]"),
    "xz": (1, "z", "x", "fz", "fx", "z [mm]", "x [mm]"),
    "yz": (0, "z", "y", "fz", "fy", "z [mm]", "y [mm]"),
}

# 旧 axis 语义 (固定/切割方向) -> 显示平面 (互为反义)
_AXIS_TO_PLANE = {"x": "yz", "y": "xz", "z": "xy"}


def resolve_3d_plane(field: FieldMap3D, plane: str = "auto",
                     mode: str = "plane") -> str:
    """归一化 plane; 'auto' 按场景选:
      * mode='plane' -> 固定网格点最少的方向 (显示平面最大, 避免薄轴色带);
      * mode='depth' -> 固定网格点最多的方向 (3D 视图体现纵深)。
    """
    if plane in _PLANE_SPEC:
        return plane
    if plane != "auto":
        raise ValueError("plane 必须为 'xy'/'xz'/'yz' 或 'auto': %r" % (plane,))
    counts = [len(field.x), len(field.y), len(field.z)]
    fixed = "xyz"[int(np.argmin(counts) if mode == "plane"
                      else np.argmax(counts))]
    return _AXIS_TO_PLANE[fixed]


def plane_fixed_axis(plane: str) -> str:
    """平面 -> 固定轴 ('xy'->'z', 'xz'->'y', 'yz'->'x'), 供交互组件复用。"""
    return "xyz"[_PLANE_SPEC[plane][0]]


def _slice_indices(n: int, n_slices: int):
    """等距取 n_slices 个剖面下标 (含首末)。"""
    if n_slices < 1:
        raise ValueError("n_slices 必须 >= 1")
    return [int(round(k * (n - 1) / max(n_slices - 1, 1)))
            for k in range(n_slices)]


def _slice_positions(fixed, n_slices: int, position=None, index=None):
    """确定剖面下标列表: position [m] / index 指定单层, 否则等距多层。"""
    if position is not None:
        return [int(np.argmin(np.abs(np.asarray(fixed) - position)))]
    if index is not None:
        i = int(index)
        if not 0 <= i < len(fixed):
            raise ValueError("index 越界: %d (范围 0..%d)"
                             % (i, len(fixed) - 1))
        return [i]
    return _slice_indices(len(fixed), n_slices)


def _plane(field: FieldMap3D, plane: str, i: int, values=None):
    """plane 方向第 i 层剖面 (i 沿固定轴):
    (px, py, v2d, u, v, xlab, ylab) — px/py 为剖面两轴网格 [m],
    v2d/u/v 均按 (len(py), len(px)) 布局 (与 pcolormesh/quiver 一致),
    u/v 为沿 px/py 方向的场分量。
    """
    if values is None:
        values = field.magnitude
    fixed_idx, pxname, pyname, uname, vname, xlab, ylab = _PLANE_SPEC[plane]
    px = getattr(field, pxname)
    py = getattr(field, pyname)
    u = getattr(field, uname)
    v = getattr(field, vname)
    if fixed_idx == 0:        # 固定 x -> yz 平面
        v2d, u2, v2 = values[i], u[i], v[i]
    elif fixed_idx == 1:      # 固定 y -> xz 平面
        v2d, u2, v2 = values[:, i], u[:, i], v[:, i]
    else:                     # 固定 z -> xy 平面
        v2d, u2, v2 = values[:, :, i], u[:, :, i], v[:, :, i]
    if fixed_idx == 2:        # xy: (px,py)=(nx,ny) -> 需 (ny,nx)
        v2d, u2, v2 = v2d.T, u2.T, v2.T
    return px, py, v2d, u2, v2, xlab, ylab


def _magnitude_label(field: FieldMap3D) -> str:
    return "|%s| [%s]" % (field.quantity, field.unit) if field.unit \
        else "|%s|" % field.quantity


def _box_aspect(field: FieldMap3D, aspect) -> tuple:
    """3D 盒子的显示比例 (返回 (rx, ry, rz))。

    'physical': 按物理尺寸 (mm 跨度), 几何真实但薄轴会被压扁;
    'equal'   : 三轴等长, 纯可读性;
    'grid'    : 按三轴网格点数比例;
    'auto'    : 物理比例 + 可读性下限 (最短轴不低于最长轴的 25%);
    (rx,ry,rz): 手动比例。
    """
    spans = np.array([np.ptp(field.x), np.ptp(field.y),
                      np.ptp(field.z)], dtype=float)
    if aspect == "physical":
        return tuple(spans)
    if aspect == "equal":
        return (1.0, 1.0, 1.0)
    if aspect == "grid":
        return (float(len(field.x)), float(len(field.y)),
                float(len(field.z)))
    if aspect == "auto":
        s = spans / max(spans.max(), 1e-30)
        s = np.maximum(s, 0.25)
        return tuple(s)
    if isinstance(aspect, (tuple, list)) and len(aspect) == 3:
        return tuple(float(a) for a in aspect)
    raise ValueError(
        "aspect 必须为 'physical'/'equal'/'grid'/'auto' 或 (rx, ry, rz): %r"
        % (aspect,))


def plot_3d_field_slices(field: FieldMap3D, plane="auto", component=None,
                         kind="heatmap", n_slices=3, position=None, index=None,
                         n_levels=12, max_arrows=400, figsize=(13, 4),
                         title=None):
    """场图剖面: 多层并排, 或单层 (position [m] / index 网格序号指定).

    component: None = |F| 幅值; 'x'/'y'/'z' = 带符号单分量。
    kind: 'heatmap' = 色块 (component=None 时叠加面内箭头 quiver) |
          'contour' = 填充等值线。
    plane='xy'/'xz'/'yz' 选择显示平面 (xz/yz 中 z 在横轴);
    'auto' 固定网格点最少的方向 (显示平面最大, 避免薄轴色带)。
    都不传时按 n_slices 等距多层。箭头确定性抽稀到 <= max_arrows 个,
    全部面板共用同一缩放; 面内分量全零时 (如 3D_Dipole 的 By 场在
    x-z 平面) 只显示底色, 这是物理事实而非错误。
    """
    plane = resolve_3d_plane(field, plane, mode="plane")
    fixed_axis = plane_fixed_axis(plane)
    fixed = getattr(field, fixed_axis)
    idxs = _slice_positions(fixed, n_slices, position, index)
    signed = component is not None
    data = field.component(component) if signed else field.magnitude
    vmax = float(np.max(np.abs(data)))
    if vmax == 0:
        warnings.warn(
            "3D 场图 %s 的%s全为零: 剖面将是空图"
            % (field.source, "分量 %s " % component if signed else ""),
            UserWarning, stacklevel=2)
        vmax = 1.0
    # 箭头统一缩放 (仅 heatmap + 幅值): 最长箭头约为面板尺度 20%
    max_arrow = 0.0
    if kind == "heatmap" and not signed:
        for i in idxs:
            _, _, _, u, v, _, _ = _plane(field, plane, i)
            m = np.sqrt(u ** 2 + v ** 2)
            if m.size:
                max_arrow = max(max_arrow, float(m.max()))
    # 面板尺度取显示平面两轴中较大者 (非全局三轴): 腔场等 z 远长于 x/y
    # 的场图若按全局跨度缩放, 箭头会被放大到溢出 xy 面板。
    _, pxname, pyname, _, _, _, _ = _PLANE_SPEC[plane]
    span_mm = max(np.ptp(getattr(field, pxname)),
                  np.ptp(getattr(field, pyname))) * 1e3
    scale = max_arrow / max(0.2 * span_mm, 1e-30) if max_arrow > 0 else None
    step_max = int(np.sqrt(max_arrows))

    fig, axs = plt.subplots(1, len(idxs), figsize=figsize, squeeze=False)
    axs = axs[0]
    if signed:
        label = ("$%s_%s$ [%s]" % (field.quantity, component, field.unit)
                 if field.unit else "$%s_%s$" % (field.quantity, component))
    else:
        label = _magnitude_label(field)
    levels = np.linspace(-vmax if signed else 0, vmax, n_levels)
    for k, i in enumerate(idxs):
        px, py, v2d, u, v, xlab, ylab = _plane(field, plane, i, values=data)
        ax = axs[k]
        if kind == "contour":
            mappable = ax.contourf(px * 1e3, py * 1e3, v2d, levels=levels,
                                   cmap="viridis" if not signed else "RdBu_r")
            ax.contour(px * 1e3, py * 1e3, v2d, levels=levels,
                       colors="k", linewidths=0.4, alpha=0.35)
        else:  # heatmap
            mappable = ax.pcolormesh(
                px * 1e3, py * 1e3, v2d,
                cmap="viridis" if not signed else "RdBu_r",
                vmin=0 if not signed else -vmax, vmax=vmax, shading="auto")
            if max_arrow > 0:
                sx = slice(None, None,
                           max(1, int(np.ceil(len(px) / step_max))))
                sy = slice(None, None,
                           max(1, int(np.ceil(len(py) / step_max))))
                ax.quiver(px[sx] * 1e3, py[sy] * 1e3, u[sy, sx], v[sy, sx],
                          angles="xy", scale_units="xy", scale=scale,
                          pivot="mid", color="0.15", width=0.003)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title("%s = %.4g mm" % (fixed_axis, fixed[i] * 1e3))
    fig.colorbar(mappable, ax=axs, label=label)
    if title:
        fig.suptitle(title)
    fig.subplots_adjust(wspace=0.35, bottom=0.15, top=0.88)
    return fig


def plot_3d_field_contour_3d(field: FieldMap3D, plane="auto", n_levels=5,
                             n_planes=7, figsize=(9, 7), title=None,
                             aspect="auto"):
    """3D 等值线: |F| 沿 plane 法向的多层 offset 平面叠加 (matplotlib 原生).

    依赖约束 (核心仅 numpy/scipy/matplotlib/pandas, 无 marching cubes)
    下用 Axes3D.contour 的多层平面展示全场强度分布; plane='auto' 固定
    网格点最多的方向以体现纵深 (如 3D_Dipole 沿 z)。真等值面
    (isosurface) 需额外依赖。
    aspect: 盒子显示比例, 见 _box_aspect (默认 'auto': 物理比例 +
    可读性下限, 避免 y=6mm 对 z=450mm 这类压扁)。
    """
    plane = resolve_3d_plane(field, plane, mode="depth")
    fixed_axis = plane_fixed_axis(plane)
    fixed = getattr(field, fixed_axis)
    mag = field.magnitude
    vmax = float(np.max(mag))
    fig = plt.figure(figsize=figsize)
    ax = fig.add_subplot(111, projection="3d")
    if vmax == 0:
        warnings.warn(
            "3D 场图 %s 的数据全为零: 无法生成 3D 等值线" % field.source,
            UserWarning, stacklevel=2)
        ax.text2D(0.5, 0.5, "field is zero", ha="center",
                  transform=ax.transAxes)
        return fig
    levels = np.linspace(0.25 * vmax, 0.95 * vmax, n_levels)
    # 3D 轴用自然坐标 (x,y,z); 剖面 = 非固定两轴
    xname, yname = [c for c in "xyz" if c != fixed_axis]
    px = getattr(field, xname)
    py = getattr(field, yname)
    PX, PY = np.meshgrid(px * 1e3, py * 1e3)
    for off in np.linspace(fixed[0], fixed[-1], n_planes):
        i = int(np.argmin(np.abs(fixed - off)))
        if fixed_axis == "z":
            v2d = mag[:, :, i]
        elif fixed_axis == "y":
            v2d = mag[:, i, :]
        else:
            v2d = mag[i, :, :]
        ax.contour(PX, PY, v2d.T, levels=levels, zdir=fixed_axis,
                   offset=off * 1e3, cmap="viridis", linewidths=1.0)
    # 3D 轴物理轴恒为 x/y/z, 标签必须写死, 不能用非固定轴名 (bug 修复)
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_zlabel("z [mm]")
    ax.set_box_aspect(_box_aspect(field, aspect))
    mappable = plt.cm.ScalarMappable(
        norm=Normalize(levels.min(), levels.max()), cmap="viridis")
    fig.colorbar(mappable, ax=ax, label=_magnitude_label(field))
    ax.set_title(title or "|%s| 3D contours along %s"
                 % (field.quantity, fixed_axis))
    return fig


def plot_3d_field_quiver(
    field: Union[FieldMap3D, str, Path],
    plane: str = "xy",
    index=None,
    position=None,
    figsize=(7, 6),
    title=None,
    color_by=None,
    max_arrows=2500,
) -> plt.Figure:
    """3D 场图 2D quiver 剖面 (纯矢量箭头图, 无底色).

    选定显示平面 plane ('xy'/'xz'/'yz'), 显示面内 2D 箭头图:
      plane='xy' (固定 z), 箭头 (U,V) = (fx, fy), 即 (Ex,Ey)/(Bx,By);
      plane='xz' (固定 y, z 在横轴), (U,V) = (fz, fx);
      plane='yz' (固定 x, z 在横轴), (U,V) = (fz, fy)。
    index (网格序号) 或 position [m] 指定单个切片位置; 都不传默认中间层。
    纯 2D quiver: angles='xy', scale_units='xy' 保证箭头方向/长度真实
    (参考 matplotlib 标准 quiver 用法, 单面板, 无 |F| 底色)。
    field 可为 FieldMap3D 或主名/分量文件路径 (自动读取, 同
    plot_3d_field_map)。
    color_by: None = 单色箭头 (C0); 'magnitude' = 按面内模长
    sqrt(U^2+V^2) 着色并加色条。
    max_arrows: 箭头确定性抽稀上限 (等步长, 超大网格防卡顿);
    默认 2500, 常见场图剖面 (11x11 / 40x46) 均不抽稀。
    """
    from ..io.field_map import FieldMap3D
    if isinstance(field, (str, Path)):
        field = read_3d_field_map_components(field)
    if not isinstance(field, FieldMap3D):
        raise TypeError("plot_3d_field_quiver 需要 FieldMap3D")
    if plane not in _PLANE_SPEC:
        raise ValueError("plane 必须为 'xy'/'xz'/'yz': %r" % (plane,))
    fixed_axis = plane_fixed_axis(plane)
    fixed = getattr(field, fixed_axis)
    if index is None and position is None:
        index = len(fixed) // 2          # 默认中间层
    i = _slice_positions(fixed, 1, position, index)[0]
    px, py, _, u, v, xlab, ylab = _plane(field, plane, i)
    if u.size > max_arrows:
        step = int(np.ceil(np.sqrt(u.size / max_arrows)))
        sx = slice(None, None, step)
        sy = slice(None, None, step)
        px, py, u, v = px[sx], py[sy], u[sy, sx], v[sy, sx]

    fig, ax = plt.subplots(figsize=figsize)
    inplane = np.sqrt(u ** 2 + v ** 2)
    if not np.any(inplane):
        # 面内分量全零 (物理事实, 如 3D_Dipole 的 B 场只在 y 方向,
        # 其 xz 平面无箭头): 跳过 quiver 以避免 matplotlib 除零,
        # 仍返回带轴/标签的图便于定位。
        warnings.warn(
            "3D 场图 %s 在 %s=%.4g mm 剖面的面内分量全为零: quiver 为空"
            % (field.source, fixed_axis, fixed[i] * 1e3),
            UserWarning, stacklevel=2)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_title(title or "%s = %.4g mm" % (fixed_axis, fixed[i] * 1e3))
        fig.tight_layout()
        return fig
    if color_by == "magnitude":
        q = ax.quiver(px * 1e3, py * 1e3, u, v, np.sqrt(u ** 2 + v ** 2),
                      angles="xy", scale_units="xy", cmap="viridis",
                      width=0.004)
        fig.colorbar(q, ax=ax, label=_magnitude_label(field))
    elif color_by is None:
        ax.quiver(px * 1e3, py * 1e3, u, v, color="C0", angles="xy",
                  scale_units="xy", width=0.004)
    else:
        raise ValueError(
            "color_by 必须为 None 或 'magnitude': %r" % (color_by,))
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    ax.set_title(title or "%s = %.4g mm" % (fixed_axis, fixed[i] * 1e3))
    fig.tight_layout()
    return fig


_VIEWS = ("slices", "contour3d", "quiver")


def plot_3d_field_map(base, view="slices", plane="auto", component=None,
                      kind="heatmap", n_slices=3, position=None, index=None,
                      n_levels=12, n_planes=7, title=None, aspect="auto"):
    """统一入口: 3D 场图剖面 / 3D 等值线 / 单层 quiver.

    Args:
        base: 场图主名 (如 3D_Dipole) 或任一分量文件 (如 3D_Dipole.by);
              自动读取全部三个分量并推导单位 (bx/by/bz -> T,
              ex/ey/ez -> V/m)。也可直接传已加载的 FieldMap3D (不再重读)。
        view: 'slices' (默认; 多层剖面, 见 component/kind) |
              'contour3d' (3D 等值线, 多层 offset 平面) |
              'quiver' (单层纯箭头)。
        plane: 'xy'/'xz'/'yz' 显示平面 (xz/yz 中 z 在横轴);
               'auto' 剖面类固定最薄方向, 3D 视图固定最密方向。
        component: None = |F| 幅值; 'x'/'y'/'z' = 带符号单分量 (仅 slices)。
        kind: 'heatmap' (默认, |F| 色块 + 面内箭头) | 'contour' (填充等值线)。
        position [m] / index: 指定单个剖面 (不传则 n_slices 等距多层)。
        n_slices / n_levels / n_planes: 剖面数 / 等值层级数 / 3D 平面数。
        aspect: 仅 contour3d; 'auto'/'physical'/'equal'/'grid' 或 (rx,ry,rz)。
    """
    if isinstance(base, FieldMap3D):
        field = base
    else:
        field = read_3d_field_map_components(base)
    if view == "slices":
        return plot_3d_field_slices(field, plane=plane, component=component,
                                    kind=kind, n_slices=n_slices,
                                    position=position, index=index,
                                    n_levels=n_levels, title=title)
    if view == "contour3d":
        return plot_3d_field_contour_3d(field, plane=plane, n_levels=n_levels,
                                        n_planes=n_planes, title=title,
                                        aspect=aspect)
    if view == "quiver":
        return plot_3d_field_quiver(field, plane=plane, position=position,
                                    index=index, title=title)
    raise ValueError("未知 view %r (可用: %s)" % (view, ", ".join(_VIEWS)))

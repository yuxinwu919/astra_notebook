"""核心发射度 (手册 4.13.5, postpro 5.6.1 项 6).

ASTRA 方法: 各粒子对 rms 发射度的贡献用**单粒子归一化相空间振幅**
    J_i = (gamma_T u_i^2 + 2 alpha u_i u'_i + beta u'_i^2) / 2
表示, 升序排序后按粒子百分比截取核心。核心发射度为核心粒子的
**平均振幅** (等价于全束团 rms 发射度减去被排除粒子的平均贡献):
    eps_core(f) = (1/N) * sum(J_i | 前 f*N 个)
由于 <J> = eps (用束团 Twiss 代入), f=1 时自动退化为标准 rms 发射度。

口径说明 (2026-08-19 R4 真跑算法破解, golden=Manual_Example/Cemit):
  * 旧实现"核心子集重算 rms 发射度"与 ASTRA 偏差 +5%/+10.5%/+25%
    (95/90/80); 均值振幅口径在 z=1.5 golden 位置三平面三分数
    max|dev| = 0.024%, z 平面全线 499 位置 max 0.001%, 横向在线圈
    区外 < 0.1%, 线圈区内 max 0.56% (ASCII 相位 dump 5 位有效数字
    舍入 + 正则动量重建噪声, 边界粒子翻转效应)。
  * 横向前置因子: J 用归一化振幅 (x/y 乘 beta*gamma = p_ref/mc),
    单位 m.rad (pi mm mrad 数值 = x1e6), 与 Cemit 文件列口径一致;
    纵向 J 为 eV.m (u = z - <z>, u' = E_kin - <E_kin>, keV.mm 数值
    上 = eV.m)。
  * 核心子集大小 k = round(f*N); N=500 时 f*N 恒为整数,
    round/floor/ceil 不可区分 (本数据集无法裁决, 保留 round)。

坐标口径与统计模块一致:
  * 横向 (x/y): u = 坐标, u' = 正则散角 (手册 4.13.1) / p_ref;
  * 纵向 (z):   u = z - <z>, u' = E_kin - <E_kin> (对应 Zemit 口径)。

加权用 |q| (混合电荷符号安全, AGENTS.md 规则 7)。

参考: ASTRA Manual V3.2 4.13.5; P. Emma 'Some basic features of the
beam emittance', PRST-AB 6, 034202 (2003) 的单粒子振幅不变量。
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from ..constants import M_E_C2_EV, kinetic_energy_from_momentum_vector
from ..distribution import Distribution
from .emittance import (
    canonical_divergence,
    canonical_signs,
    compute_geometric_emittance,
    compute_twiss_parameters,
)


def _plane_coords(dist: Distribution, plane: str, bz_on_axis_T: float,
                  ref_momentum_eVc: float, weights: Optional[np.ndarray]):
    """返回 (u, up, w): 居中坐标与散角 (无单位), |q| 权重。

    横向: up = 正则散角 / p_ref (无单位, 与 Xemit 一致)。
    纵向: up = E_kin - <E_kin> [eV] (与 Zemit / emit_z 口径一致)。
    """
    m = dist.active
    x = dist.x[m]
    y = dist.y[m]
    pz = dist.pz[m]
    w: Optional[np.ndarray] = (np.abs(weights[m])
                               if weights is not None else None)

    def _avg(a: np.ndarray) -> float:
        return float(np.average(a, weights=w) if w is not None else np.mean(a))

    if plane == "x":
        s_can = canonical_signs(dist)[m]   # 种类感知符号 (2026-08 F4)
        ptx = canonical_divergence(dist.px[m], y, bz_on_axis_T, s_can)
        u = x - _avg(x)
        up = (ptx - _avg(ptx)) / ref_momentum_eVc
    elif plane == "y":
        s_can = canonical_signs(dist)[m]
        pty = canonical_divergence(dist.py[m], x, bz_on_axis_T, -s_can)
        u = y - _avg(y)
        up = (pty - _avg(pty)) / ref_momentum_eVc
    elif plane == "z":
        # 全动量动能 (2026-08 审计 P2-2)
        e = kinetic_energy_from_momentum_vector(dist.px[m], dist.py[m], pz)
        u = dist.z[m] - _avg(dist.z[m])
        up = e - _avg(e)
    else:
        raise ValueError("plane 必须为 'x'/'y'/'z': %r" % (plane,))
    return np.asarray(u, dtype=float), np.asarray(up, dtype=float), w


def single_particle_amplitudes(
    dist: Distribution,
    plane: str = "x",
    bz_on_axis_T: float = 0.0,
    ref_momentum_eVc: Optional[float] = None,
) -> np.ndarray:
    """单粒子发射度不变量 J_i (手册 4.13.5), 按 active 粒子排列.

    横向 (x/y): 归一化振幅 [m.rad] (pi mm mrad 数值 = x1e6, 与 Cemit
    文件列同口径, J x beta*gamma = p_ref/mc);
    纵向 (z): [eV.m] (u=z, u'=E, keV.mm 数值 = x1)。
    """
    if not np.any(dist.active):
        raise ValueError("no active particles (status>1)")
    if ref_momentum_eVc is None:
        ref_momentum_eVc = dist.ref_momentum_or_mean()
    if ref_momentum_eVc <= 0:
        raise ValueError("reference momentum is zero/negative")
    u, up, w = _plane_coords(dist, plane, bz_on_axis_T, ref_momentum_eVc,
                             dist.charge)
    eps = compute_geometric_emittance(u, up, w)
    if eps <= 0:
        return np.zeros(len(u))
    beta, alpha, gamma_t = compute_twiss_parameters(u, up, w)
    j = 0.5 * (gamma_t * u ** 2 + 2.0 * alpha * u * up + beta * up ** 2)
    if plane in ("x", "y"):
        # ASTRA Cemit 横向列是归一化发射度 (R4 实测): 振幅乘 beta*gamma.
        j = j * (ref_momentum_eVc / M_E_C2_EV)
    return j


def compute_core_emittance_by_fraction(
    dist: Distribution,
    plane: str = "x",
    fractions=(1.0, 0.95, 0.90, 0.80),
    bz_on_axis_T: float = 0.0,
    ref_momentum_eVc: Optional[float] = None,
) -> dict:
    """核心发射度 vs 粒子百分比 (手册 4.13.5 口径).

    ASTRA 算法 (R4 真跑破译, 见模块 docstring): 按单粒子振幅 J_i 升序
    排序后, 核心发射度 = 前 f*N 个粒子的平均振幅 sum(J)/N — 不是核心
    子集重算 rms 发射度 (旧口径偏差 +5/+10.5/+25%)。

    Returns:
        {fraction: eps} — 横向 eps_n [m.rad] (pi mm mrad 数值 = x1e6),
        纵向 eps_z [eV.m] (keV.mm 数值 = x1)。与 read_cemit_file 一致。
    """
    if not np.any(dist.active):
        raise ValueError("no active particles (status>1)")
    if ref_momentum_eVc is None:
        ref_momentum_eVc = dist.ref_momentum_or_mean()
    if ref_momentum_eVc <= 0:
        raise ValueError("reference momentum is zero/negative")

    j = single_particle_amplitudes(dist, plane, bz_on_axis_T, ref_momentum_eVc)
    order = np.argsort(j)

    out = {}
    n = len(order)
    for f in fractions:
        if f <= 0:
            out[float(f)] = 0.0
            continue
        k = max(int(round(f * n)), 1)
        sel = order[:k]
        out[float(f)] = float(np.sum(j[sel]) / n)
    return out


def compute_core_emittance_curves(
    dist: Distribution,
    planes=("x", "y", "z"),
    fractions=(0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0),
    bz_on_axis_T: float = 0.0,
) -> dict:
    """三平面核心发射度曲线, 供绘图. {plane: {fraction: eps}}."""
    return {plane: compute_core_emittance_by_fraction(
        dist, plane=plane, fractions=fractions, bz_on_axis_T=bz_on_axis_T)
        for plane in planes}

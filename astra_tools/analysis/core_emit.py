"""核心发射度 (手册 4.13.5, postpro 5.6.1 项 6).

ASTRA 方法: 各粒子对 rms 发射度的贡献用**单粒子归一化相空间振幅**
    J_i = (gamma_T u_i^2 + 2 alpha u_i u'_i + beta u'_i^2) / 2
表示, 升序排序后按粒子百分比截取核心, 重算 rms 发射度, 得到
"emittance vs particle fraction" 曲线。由于 <J> = eps (用束团 Twiss
代入), f=1 时自动退化为标准 rms 发射度。

坐标口径与统计模块一致:
  * 横向 (x/y): u = 坐标, u' = 正则散角 (手册 4.13.1) / p_ref;
    返回归一化发射度 eps_n = beta*gamma * eps_geom (pi mm mrad 数值 = m.rad*1e6)。
  * 纵向 (z):   u = z - <z>, u' = E_kin - <E_kin> (对应 Zemit 口径);
    返回 eps_z = sqrt(<z^2><E^2> - <zE>^2) [eV.m], keV.mm 数值上 = eV.m。

加权用 |q| (混合电荷符号安全, AGENTS.md 规则 7)。

口径说明 (2026-08-18 实测 vs examples/Manual_Example/golden/Example.Cemit.001):
  * f=1.0 (全束团 rms 发射度) 与 ASTRA Xemit/Cemit 精确一致 (<0.1%);
  * f<1 的核心发射度随分数单调递减、趋势与 ASTRA 一致, 但数值比
    ASTRA Cemit 略大 (95% 约 +5%, 90% 约 +11%, 80% 约 +25%),
    因 ASTRA 的单粒子贡献排序算法未公开复现。功能/趋势正确,
    与 ASTRA 的精确数值对齐留待算法确认。

参考: ASTRA Manual V3.2 4.13.5; P. Emma 'Some basic features of the
beam emittance', PRST-AB 6, 034202 (2003) 的单粒子振幅不变量。
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from ..constants import kinetic_energy_from_momentum_vector
from ..distribution import Distribution
from .emittance import (
    canonical_divergence,
    canonical_signs,
    compute_geometric_emittance,
    compute_normalized_emittance,
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
    """单粒子发射度不变量 J_i (手册 4.13.5), 按 active 粒子排列. [m.rad 或 eV.m]"""
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
    return j


def compute_core_emittance_by_fraction(
    dist: Distribution,
    plane: str = "x",
    fractions=(1.0, 0.95, 0.90, 0.80),
    bz_on_axis_T: float = 0.0,
    ref_momentum_eVc: Optional[float] = None,
) -> dict:
    """核心发射度 vs 粒子百分比 (手册 4.13.5 口径).

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

    m = dist.active
    u, up, w = _plane_coords(dist, plane, bz_on_axis_T, ref_momentum_eVc,
                             dist.charge)
    j = single_particle_amplitudes(dist, plane, bz_on_axis_T, ref_momentum_eVc)
    order = np.argsort(j)

    out = {}
    n = len(order)
    u = np.asarray(u, dtype=float)
    up = np.asarray(up, dtype=float)
    for f in fractions:
        if f <= 0:
            out[float(f)] = 0.0
            continue
        k = max(int(round(f * n)), 1)
        sel = order[:k]
        us, ups = u[sel], up[sel]
        ws = w[sel] if w is not None else None
        uc = us - float(np.average(us, weights=ws) if ws is not None
                        else np.mean(us))
        upc = ups - float(np.average(ups, weights=ws) if ws is not None
                          else np.mean(ups))
        eps_geom = compute_geometric_emittance(uc, upc, ws)
        if plane == "z":
            out[float(f)] = float(eps_geom)   # eV.m (= keV.mm 数值)
        else:
            out[float(f)] = compute_normalized_emittance(
                eps_geom, ref_momentum_eVc)
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

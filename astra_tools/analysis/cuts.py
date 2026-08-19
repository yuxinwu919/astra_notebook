"""Phase space cuts and rotations (postpro 5.6.4).

Apply transverse/longitudinal/energy windows or a radial aperture to a
Distribution, or rotate the x-y plane; returns a new Distribution ready
for further tracking or export.
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from ..constants import kinetic_energy_from_momentum_vector
from ..distribution import Distribution


def cut_distribution(
    dist: Distribution,
    x_range=None,
    y_range=None,
    z_range=None,
    e_range=None,
    r_aperture=None,
):
    """Cut a distribution by windows (active particles only; cut
    particles are relabelled status -31, 'discarded by user', per the
    ASTRA Manual Table 2).

    Args:
        dist: input distribution.
        x_range / y_range: (min, max) [m] transverse windows.
        z_range: (min, max) [m] longitudinal window (absolute z).
        e_range: (min, max) [eV] kinetic-energy window.
        r_aperture: radial aperture radius [m], sqrt(x^2+y^2) < r.

    Returns:
        (dist_cut, mask) - mask marks the removed particles.
    """
    cut = np.zeros(dist.n_particle, dtype=bool)
    if x_range is not None:
        cut |= (dist.x < x_range[0]) | (dist.x > x_range[1])
    if y_range is not None:
        cut |= (dist.y < y_range[0]) | (dist.y > y_range[1])
    if z_range is not None:
        cut |= (dist.z < z_range[0]) | (dist.z > z_range[1])
    if e_range is not None:
        e = kinetic_energy_from_momentum_vector(dist.px, dist.py, dist.pz)
        cut |= (e < e_range[0]) | (e > e_range[1])
    if r_aperture is not None:
        r = np.sqrt(dist.x**2 + dist.y**2)
        cut |= r > r_aperture

    cut &= dist.active
    status = dist.status.copy()
    status[cut] = -31
    out = Distribution(
        x=dist.x.copy(), y=dist.y.copy(), z=dist.z.copy(),
        px=dist.px.copy(), py=dist.py.copy(), pz=dist.pz.copy(),
        clock=dist.clock.copy(), charge=dist.charge.copy(),
        status=status,
        index=None if dist.index is None else dist.index.copy(),
        ref_time_ns=dist.ref_time_ns, ref_momentum_eVc=dist.ref_momentum_eVc,
        total_charge_nC=dist.total_charge_nC,
        ref_x_m=dist.ref_x_m, ref_y_m=dist.ref_y_m, ref_z_m=dist.ref_z_m,
        source=dist.source + " (cut)", format=dist.format, attrs=dict(dist.attrs),
    )
    return out, cut


def rotate_phase_space(dist: Distribution, angle_deg: float):
    """Rotate the x-y plane (coordinates and momenta together)."""
    th = np.radians(angle_deg)
    c, s = np.cos(th), np.sin(th)
    out = Distribution(
        x=c * dist.x + s * dist.y,
        y=-s * dist.x + c * dist.y,
        z=dist.z.copy(),
        px=c * dist.px + s * dist.py,
        py=-s * dist.px + c * dist.py,
        pz=dist.pz.copy(), clock=dist.clock.copy(), charge=dist.charge.copy(),
        status=dist.status.copy(),
        index=None if dist.index is None else dist.index.copy(),
        ref_time_ns=dist.ref_time_ns, ref_momentum_eVc=dist.ref_momentum_eVc,
        total_charge_nC=dist.total_charge_nC,
        ref_x_m=dist.ref_x_m, ref_y_m=dist.ref_y_m, ref_z_m=dist.ref_z_m,
        source=dist.source + " (rotated %.1f deg)" % angle_deg,
        format=dist.format, attrs=dict(dist.attrs),
    )
    return out


def modify_correlated_energy_spread(dist: Distribution, factor: float):
    """改变关联能散 (postpro 5.6.3 项 12).

    把 pz 分解为"非相关部分 + 随 z 线性相关部分":
        pz = pz_uncorr + (a + b*z)
    乘以 factor 缩放关联部分并重新合成:
        pz' = pz_uncorr + factor * (a + b*z)
    factor=1 不变, 0 完全去相关, >1 增强关联。
    返回新 Distribution (active 粒子参与拟合, 全体粒子变换)。
    """
    m = dist.active
    z = dist.z
    pz = dist.pz
    if m.sum() < 3:
        raise ValueError("too few active particles for correlation fit")
    b, a = np.polyfit(z[m], pz[m], 1)
    corr = a + b * z
    pz_new = (pz - corr) + factor * corr
    return Distribution(
        x=dist.x.copy(), y=dist.y.copy(), z=z.copy(),
        px=dist.px.copy(), py=dist.py.copy(), pz=pz_new,
        clock=dist.clock.copy(), charge=dist.charge.copy(),
        status=dist.status.copy(),
        index=None if dist.index is None else dist.index.copy(),
        ref_time_ns=dist.ref_time_ns, ref_momentum_eVc=dist.ref_momentum_eVc,
        total_charge_nC=dist.total_charge_nC,
        ref_x_m=dist.ref_x_m, ref_y_m=dist.ref_y_m, ref_z_m=dist.ref_z_m,
        source=dist.source + " (corr. E spread x%.2f)" % factor,
        format=dist.format, attrs=dict(dist.attrs),
    )


def optimized_cut_center(values: np.ndarray, width: float,
                         weights: Optional[np.ndarray] = None) -> float:
    """优化切割中心 (手册 5.6.4): 找宽度为 width 的对称窗口中心 c,
    使窗口 [c-width/2, c+width/2] 内 (|q| 加权) 计数最大."""
    v = np.asarray(values, dtype=float)
    w = np.abs(np.asarray(weights, dtype=float)) if weights is not None else None
    if len(v) == 0 or width <= 0:
        raise ValueError("empty values or non-positive width")
    order = np.argsort(v)
    vs = v[order]
    ws = w[order] if w is not None else np.ones_like(vs)
    # 前缀和, 滑动窗口最大计数
    cum = np.concatenate([[0.0], np.cumsum(ws)])
    best, best_c = -1.0, float(np.mean(vs))
    lo = 0
    for hi in range(len(vs)):
        while vs[lo] < vs[hi] - width:
            lo += 1
        total = cum[hi + 1] - cum[lo]
        if total > best:
            best = total
            best_c = 0.5 * (vs[lo] + vs[hi])
    return float(best_c)


def optimized_cut(dist: Distribution, width: float, param: str = "z"):
    """优化切割 (手册 5.6.4): 给定区间参数 width (如束长), 在 active
    粒子中找到使存活 (|q| 加权) 粒子数最大的对称窗口, 返回切割分布.

    param: 'x'/'y'/'z' (SI [m]) 或 'E' (动能 [eV], width 用 [eV])。
    """
    m = dist.active
    if param == "x":
        values = dist.x
    elif param == "y":
        values = dist.y
    elif param == "z":
        values = dist.z
    elif param == "E":
        values = kinetic_energy_from_momentum_vector(dist.px, dist.py, dist.pz)
    else:
        raise ValueError("param 必须为 'x'/'y'/'z'/'E': %r" % (param,))
    c = optimized_cut_center(np.asarray(values[m], dtype=float), width,
                             dist.charge[m])
    lo, hi = c - 0.5 * width, c + 0.5 * width
    return cut_distribution(dist, **{param + "_range": (lo, hi)})

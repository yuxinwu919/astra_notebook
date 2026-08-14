"""核心电荷分数分析 (core charge-fraction curves).

postpro 菜单的"核心束长/发射度 vs 电荷分数"曲线: 把束团按纵向 z 排序,
取中间 charge_fraction 的电荷作为"核心", 计算其束长 / 尺寸 / 归一化
发射度。charge_fraction = 1.0 即整个束团 (与 compute_statistics 一致,
可作自检点)。

约定与 statistics.py 相同: 活粒子 (status>1)、群体矩 (ddof=0)、
归一化发射度用参考动量 p_ref 与正则动量 (bz_on_axis_T)。
"""

from __future__ import annotations

from typing import Optional

import numpy as np

from ..constants import C_LIGHT, gamma_from_momentum
from ..distribution import Distribution
from .emittance import compute_geometric_emittance


def compute_core_fraction_curves(
    dist: Distribution,
    fractions=(0.1, 0.25, 0.5, 0.75, 0.9, 1.0),
    ref_momentum_eVc: Optional[float] = None,
    bz_on_axis_T: float = 0.0,
) -> dict:
    """按纵向电荷分数取核心并计算参数曲线。

    Args:
        dist: 粒子分布 (使用活粒子)。
        fractions: 核心电荷分数序列 (0..1)。
        ref_momentum_eVc: 归一化发射度参考动量; 默认取分布头参考动量,
            否则活粒子平均 pz。
        bz_on_axis_T: 正则动量的螺线管轴上场 [T]。

    Returns:
        dict: fractions, sig_z [m], sig_x [m], sig_y [m],
              emit_xn [m.rad], emit_yn [m.rad], n_particles。
    """
    d = dist.filter_active()
    if d.n_active < 6:
        raise ValueError("too few active particles for core analysis")
    if ref_momentum_eVc is None:
        ref_momentum_eVc = (
            d.ref_momentum_eVc if d.ref_momentum_eVc != 0
            else float(np.mean(d.pz))
        )
    if ref_momentum_eVc <= 0:
        raise ValueError(
            "reference momentum is zero/negative (beam at rest?); "
            "core emittances are undefined")

    order = np.argsort(d.z)
    q_abs = np.abs(d.charge[order])
    cum = np.cumsum(q_abs)
    q_tot = float(cum[-1])

    bg = float(np.sqrt(max(gamma_from_momentum(ref_momentum_eVc) ** 2 - 1.0, 0.0)))
    n_core = np.zeros(len(fractions), dtype=int)
    sig_z = np.zeros(len(fractions))
    sig_x = np.zeros(len(fractions))
    sig_y = np.zeros(len(fractions))
    emit_xn = np.zeros(len(fractions))
    emit_yn = np.zeros(len(fractions))

    for i, f in enumerate(fractions):
        if f >= 1.0:
            idx = order  # 全束团 (仍按 z 排序, 不影响统计)
        else:
            lo = (1.0 - f) / 2.0 * q_tot
            hi = q_tot - lo
            idx = order[(cum >= lo) & (cum <= hi)]
        if len(idx) < 3:
            continue
        n_core[i] = len(idx)
        xi = d.x[idx]
        yi = d.y[idx]
        ptx = d.px[idx] + 0.5 * C_LIGHT * bz_on_axis_T * yi
        pty = d.py[idx] - 0.5 * C_LIGHT * bz_on_axis_T * xi
        sig_z[i] = float(np.std(d.z[idx]))
        sig_x[i] = float(np.std(xi))
        sig_y[i] = float(np.std(yi))
        xp = (ptx - np.mean(ptx)) / ref_momentum_eVc
        yp = (pty - np.mean(pty)) / ref_momentum_eVc
        # 群体矩 (无加权), 与 compute_statistics 默认约定一致
        emit_xn[i] = bg * compute_geometric_emittance(xi - np.mean(xi), xp)
        emit_yn[i] = bg * compute_geometric_emittance(yi - np.mean(yi), yp)

    return {
        "fractions": np.asarray(fractions, dtype=float),
        "n_particles": n_core,
        "sig_z": sig_z, "sig_x": sig_x, "sig_y": sig_y,
        "emit_xn": emit_xn, "emit_yn": emit_yn,
    }

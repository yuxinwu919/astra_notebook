"""时间坐标 (postpro 5.6 菜单 1 项 3/5: 纵向相空间·时间, 三视图 vs 时间).

手册 5.6: 以时间为坐标的图显示——
  * 未发射粒子 (status -1..-6): 发射时间 (clock);
  * 已发射粒子 (status >= 0): 由相对纵向位置和束团平均纵向速度计算
        t = (z - <z>) / (beta_bar * c);
  * 混合分布 (既有未发射又有已发射粒子) 给出警告。

所有返回为 SI (秒)。
"""

from __future__ import annotations

import warnings

import numpy as np

from ..constants import C_LIGHT, beta_from_gamma, gamma_from_momentum
from ..distribution import Distribution


def bunch_time(dist: Distribution, warn_mixed: bool = True) -> np.ndarray:
    """返回每个粒子的时间坐标 [s] (手册 5.6 规则).

    Args:
        dist: 粒子分布。
        warn_mixed: 混合分布 (未发射 + 已发射同时存在) 是否发警告。

    Returns:
        t [s] (与 dist 等长)。未发射粒子用 clock, 其余用
        (z - <z>_active) / (beta_bar * c)。
    """
    n = dist.n_particle
    t = np.zeros(n, dtype=float)

    active = dist.active
    not_started = dist.not_started
    tracked = dist.tracked  # status >= 0

    # 束团平均纵向速度 (active 粒子, 手册 5.6 "average longitudinal
    # velocity of the bunch")
    beta_bar = 1.0
    if dist.n_active:
        pz_avr = float(np.mean(dist.pz[active]))
        if pz_avr > 0:
            beta_bar = beta_from_gamma(gamma_from_momentum(pz_avr))
    z_avr = float(np.mean(dist.z[active])) if dist.n_active else 0.0

    m_tracked = tracked & ~active  # passive probes also use velocity
    if np.any(tracked):
        t[tracked] = (dist.z[tracked] - z_avr) / (beta_bar * C_LIGHT)
    # passive 粒子若 z 语义不同 (探针), 仍按同一平均速度; 保留上面统一处理
    if np.any(not_started):
        t[not_started] = dist.clock[not_started]

    if warn_mixed and np.any(not_started) and np.any(m_tracked | active):
        warnings.warn(
            "混合分布: 同时存在未发射 (status -1..-6, 用 clock) 与已发射 "
            "粒子 (用束团平均速度), 时间坐标混用两种约定", UserWarning,
            stacklevel=2)
    return t


def bunch_time_ps(dist: Distribution, warn_mixed: bool = True) -> np.ndarray:
    """时间坐标 [ps] 便捷版。"""
    return bunch_time(dist, warn_mixed=warn_mixed) * 1e12

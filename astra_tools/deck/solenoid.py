"""从 ASTRA deck 解析螺线管参数并计算束团中心处轴上场 (批 3).

手册 4.13.1: 发射度应使用束团中心处螺线管轴上场构成正则动量;
此前前端从不传 bz_on_axis_T, 螺线管束团显示 trace-space 值且无告警。
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import numpy as np
from scipy.interpolate import interp1d

from ..namelist.parse import parse_namelists


def _first(v):
    """数组参数取第一个元素 (File_Bfield(1)/MaxB(1)/S_pos(1))。"""
    return v[0] if isinstance(v, (list, tuple)) else v


def solenoid_bz_at_z(
    deck_path,
    z_m: float,
    field_dir=None,
) -> Optional[float]:
    """束团中心 z 处的螺线管轴上场 [T]。

    从 deck 的 &SOLENOID 解析 File_Bfield(1)/MaxB(1)/S_pos(1),
    按场表插值: Bz(z) = interp(z - s_pos, table) * MaxB / max(|table|)。
    deck 无螺线管 (LBField!=T) 或任何解析失败时返回 None (调用方降级
    为告警, 不自动给错误数值)。
    """
    deck = Path(deck_path)
    if not deck.exists():
        return None
    try:
        blocks = parse_namelists(deck)
    except Exception:
        return None
    sol = blocks.get("SOLENOID")
    if not sol or not sol.get("LBField", False):
        return None
    try:
        fname = str(_first(sol.get("File_Bfield", ""))).strip("'").strip(chr(34))
        maxB = float(_first(sol.get("MaxB", 0.0)))
        s_pos = float(_first(sol.get("S_pos", 0.0)))
    except (TypeError, ValueError):
        return None
    if not fname or maxB <= 0:
        return None
    base = Path(field_dir) if field_dir else deck.parent
    fpath = base / fname
    if not fpath.exists():
        return None
    try:
        table = np.loadtxt(fpath, ndmin=2, encoding="utf-8")
        if table.shape[1] < 2:
            return None
    except Exception:
        return None
    zcol, bcol = table[:, 0], table[:, 1]
    bmax = float(np.max(np.abs(bcol)))
    if bmax == 0:
        return None
    val = float(interp1d(
        zcol, bcol, bounds_error=False,
        fill_value=(bcol[0], bcol[-1]))(z_m - s_pos))
    return val * maxB / bmax

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


def _as_list(v):
    """参数可能是标量或数组 (MaxB(1)/MaxB(2)...)。"""
    return list(v) if isinstance(v, (list, tuple)) else [v]


def solenoid_bz_at_z(
    deck_path,
    z_m: float,
    field_dir=None,
) -> Optional[float]:
    """束团中心 z 处的螺线管轴上场 [T] (全部螺线管叠加).

    从 deck 的 &SOLENOID 解析 File_Bfield(n)/MaxB(n)/S_pos(n),
    逐个元素按场表插值: Bz_i(z) = interp(z - s_pos_i, table_i) *
    MaxB_i / max(|table_i|), 求和 (ASTRA 追踪时所有螺线管场叠加)。
    批 6: 此前只取第一个元素, 多螺线管束线会静默用错场。
    deck 无螺线管 (LBField!=T)、任何解析失败、或某个已声明的螺线管
    (非空 File_Bfield 且 MaxB>0) 场表缺失/不可读时返回 None (调用方
    降级为告警, 不自动给错误数值、也不静默部分求和)。
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
    fnames = _as_list(sol.get("File_Bfield", ""))
    maxbs = _as_list(sol.get("MaxB", 0.0))
    sposs = _as_list(sol.get("S_pos", 0.0))
    n = max(len(fnames), len(maxbs), len(sposs))
    base = Path(field_dir) if field_dir else deck.parent
    total = 0.0
    found_any = False
    for i in range(n):
        if i >= len(fnames):
            continue
        try:
            fname = str(fnames[i]).strip("'").strip(chr(34))
            maxB = float(maxbs[min(i, len(maxbs) - 1)])
            s_pos = float(sposs[min(i, len(sposs) - 1)])
        except (TypeError, ValueError, IndexError):
            continue
        if not fname or maxB <= 0:
            continue
        fpath = base / fname
        if not fpath.exists():
            return None
        try:
            table = np.loadtxt(fpath, ndmin=2, encoding="utf-8")
        except Exception:
            return None
        if table.shape[1] < 2:
            return None
        zcol, bcol = table[:, 0], table[:, 1]
        bmax = float(np.max(np.abs(bcol)))
        if not np.isfinite(bmax):
            return None
        if bmax == 0:
            continue
        val = float(interp1d(
            zcol, bcol, bounds_error=False,
            fill_value=(bcol[0], bcol[-1]))(z_m - s_pos))
        total += val * maxB / bmax
        found_any = True
    return total if found_any else None

"""Plot_steering.par 解析 (手册 5.6: 重定向 postpro 统计/显示条件).

文件含 namelist 'Steering_parameters' (手册 5.6 / 6.1):
  * Stat(1, -100:100) / Stat(2, -100:100): 逻辑数组, 重定向显示条件。
      第 1 下标 = 2 只影响 slice/core/相空间切割图;
      第 1 下标 = 1 影响除 z-plot 外的所有图。
      某 status 对应项为 T 时才绘制该状态粒子。
  * ion_mass(1..): 自定义粒子质量 (与 NEWRUN 同义, 按电荷态归一)。
  * CP_ind_1 .. CP_ind_15: 粒子索引 -> RGB (0..1), 混合分布按粒子
      类型着色 (仅 Plot_mode=1)。RGB 为 0,0,0 (黑) 的粒子不绘制。
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np


def read_plot_steering(path) -> dict:
    """解析 Plot_steering.par, 返回:
        {"stat": {1: {flag: bool}, 2: {flag: bool}},
         "ion_mass": [...],
         "cp_ind": {idx: (r, g, b)}}
    无 Steering_parameters 块时返回空 dict。
    """
    text = Path(path).read_text(encoding="utf-8")
    m = re.search(r"&\s*Steering_parameters(.*?)/", text, re.S | re.I)
    if not m:
        return {}
    out: dict = {"stat": {}, "ion_mass": [], "cp_ind": {}}
    for line in m.group(1).splitlines():
        s = line.strip().rstrip(",")
        if not s or s.startswith("!"):
            continue
        m2 = re.match(r"Stat\s*\(\s*(\d)\s*,\s*(-?\d+)\s*\)\s*=\s*(T|F)",
                      s, re.I)
        if m2:
            idx, flag = int(m2.group(1)), int(m2.group(2))
            val = m2.group(3).upper() == "T"
            out["stat"].setdefault(idx, {})[flag] = val
            continue
        # 一行多个 Stat(...)=T/F 时逐个收集
        for m2b in re.finditer(
                r"Stat\s*\(\s*(\d)\s*,\s*(-?\d+)\s*\)\s*=\s*(T|F)",
                s, re.I):
            idx, flag = int(m2b.group(1)), int(m2b.group(2))
            val = m2b.group(3).upper() == "T"
            out["stat"].setdefault(idx, {})[flag] = val
        m3 = re.match(r"CP_ind_(\d+)\s*=\s*(.+)", s, re.I)
        if m3:
            n = int(m3.group(1))
            vals = [float(x) for x in
                    re.findall(r"[-+]?\d*\.?\d+", m3.group(2))]
            if len(vals) >= 3:
                out["cp_ind"][n] = tuple(vals[:3])
            continue
        m4 = re.match(r"ion_mass\s*\(?\s*\d*\s*\)?\s*=\s*(.+)", s, re.I)
        if m4:
            out["ion_mass"].extend(
                float(x) for x in re.findall(r"[-+]?\d*\.?\d+", m4.group(1)))
    return out


def cp_index_colors(index: np.ndarray, cp_ind: dict,
                    fallback: str = "#0077BB"):
    """混合分布按粒子索引着色 (手册 5.6 CP_ind, Plot_mode=1).

    返回 matplotlib 颜色数组; 未指定索引的粒子用 fallback;
    RGB 为黑的粒子返回 None (调用方应跳过/不绘制)。
    """
    idx = np.asarray(index)
    colors = np.array([fallback] * len(idx), dtype=object)
    for n, (r, g, b) in cp_ind.items():
        m = idx == n
        k = int(np.count_nonzero(m))
        if not k:
            continue
        if r == 0.0 and g == 0.0 and b == 0.0:
            # 黑色 = 不绘制: 标 None (调用方据此排除)
            sub = np.empty(k, dtype=object)
            sub.fill(None)
            colors[m] = sub
            continue
        sub = np.empty(k, dtype=object)
        sub.fill((r, g, b))
        colors[m] = sub
    return colors

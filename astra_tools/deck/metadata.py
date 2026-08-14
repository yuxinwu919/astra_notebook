"""ASTRA/Generator 参数元数据库 (提取自 ASTRA 手册 V3.2 第 6/7 章).

每个 namelist 一个 JSON 文件 (data/parameters/<NAME>.json), 条目:
    name    参数名 (数组参数带 () 后缀, 如 MaxE())
    type    Fortran 类型 (Real*8 / Integer / Logical / Character*150 / ... array)
    unit    单位 (来自手册; 'π mrad mm' 中的 π 表示椭圆面积语义, 数值上即 mm.mrad)
    default 默认值 (手册值; None 表示必填/无默认)
    desc    手册描述 (英文, 已清理 LaTeX)
"""

from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

_DATA_DIR = Path(__file__).resolve().parent.parent / "data" / "parameters"

# ASTRA 全部 13 个 namelist + Generator INPUT (手册第 6/7 章)
NAMELISTS = (
    "NEWRUN", "OUTPUT", "SCAN", "MODULES", "ERROR", "CHARGE",
    "APERTURE", "WAKE", "CAVITY", "SOLENOID", "QUADRUPOLE",
    "DIPOLE", "LASER", "INPUT",
)


@lru_cache(maxsize=None)
def load(namelist: str) -> dict:
    """加载一个 namelist 的完整元数据 (缓存)."""
    name = namelist.upper()
    if name not in NAMELISTS:
        raise KeyError("unknown namelist: " + name)
    path = _DATA_DIR / (name + ".json")
    with open(path, encoding="utf-8") as f:
        return json.load(f)


def params(namelist: str):
    """参数条目列表."""
    return load(namelist)["params"]


def param_info(name: str):
    """按参数名查找 (name, entry, namelist); 跨 namelist 搜索."""
    base = name.rstrip("() ")
    for nl in NAMELISTS:
        for p in params(nl):
            if p["name"].rstrip("() ") == base:
                return nl, p
    return None, None


def summary() -> str:
    """全部 namelist 的参数数概览."""
    lines = []
    for nl in NAMELISTS:
        d = load(nl)
        lines.append("  %-11s (%s) %d 参数" % (nl, d["section"], d["n_params"]))
    return "\n".join(lines)

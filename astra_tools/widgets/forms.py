"""元数据驱动的参数表单生成器 (ipywidgets).

从 astra_tools.deck.metadata 的参数元数据库 (类型/单位/默认值/描述)
自动渲染一个 namelist 的参数表单; 用户无需写代码即可填参数。
"""

from __future__ import annotations

import ipywidgets as widgets
from IPython.display import display

from ..deck.metadata import load


def _widget_for_param(entry, value=None):
    """按参数元数据类型构造控件."""
    name = entry["name"]
    ptype = entry["type"]
    desc = entry["desc"]
    unit = entry["unit"] or ""
    tip = desc[:120] + ("  [" + unit + "]" if unit else "")

    if "Logical" in ptype:
        default = False
        if value is not None:
            default = bool(value)
        elif entry["default"] and entry["default"].upper() in ("T", "TRUE"):
            default = True
        w = widgets.Checkbox(value=default, description=name, tooltip=tip)
    elif "Integer" in ptype and "array" not in ptype.lower():
        default = 0
        if value is not None:
            default = int(value)
        elif entry["default"] not in (None, ""):
            try:
                default = int(float(entry["default"]))
            except ValueError:
                default = 0
        w = widgets.IntText(value=default, description=name, tooltip=tip)
    elif "Real" in ptype and "array" not in ptype.lower():
        default = 0.0
        if value is not None:
            default = float(value)
        elif entry["default"] not in (None, ""):
            try:
                default = float(entry["default"])
            except ValueError:
                default = 0.0
        w = widgets.FloatText(value=default, description=name, tooltip=tip)
    else:
        # Character / 数组 / 其他 -> 文本框 (数组用逗号分隔)
        default = ""
        if value is not None:
            default = str(value)
        elif entry["default"] not in (None, ""):
            default = str(entry["default"])
        w = widgets.Text(value=default, description=name, tooltip=tip,
                         layout=widgets.Layout(width="70%"))
    return w


def namelist_form(namelist: str, values=None, only=None, show: bool = True):
    """为一个 namelist 生成参数表单.

    Args:
        namelist: namelist 名 (如 'INPUT', 'NEWRUN').
        values: 可选, 现有参数字典 (编辑模式).
        only: 可选, 只显示这些参数名 (常用参数模式).
        show: 是否立即 display.

    Returns:
        (widgets_dict, getter) — getter() 返回 {param: value} 字典
        (未填写的可选参数被跳过)。
    """
    meta = load(namelist)
    wmap = {}
    touched = set()
    only_lower = None
    if only is not None:
        only_lower = [o.rstrip("() ").lower() for o in only]
    for entry in meta["params"]:
        base = entry["name"].rstrip("() ")
        if only_lower is not None and base.lower() not in only_lower:
            continue
        val = None
        if values:
            for k, v in values.items():
                if k.rstrip("() ").lower() == base.lower():
                    val = v
                    break
        wmap[base] = _widget_for_param(entry, val)
        # 记录用户改动 (changed_only 模式只写用户动过的参数)
        wmap[base].observe(lambda change, b=base: touched.add(b), "value")

    box = widgets.VBox(list(wmap.values()))
    if show:
        display(widgets.HTML("<h3>&amp;%s</h3>" % namelist))
        display(box)

    def getter(changed_only: bool = False):
        out = {}
        for base, w in wmap.items():
            if changed_only and base not in touched:
                continue
            v = w.value
            # None/空串 = 未设置 (跳过); False 是语义值必须写出
            # (例如 Cathode=F 若省略会退回手册默认 T, 行为完全不同)
            if v is None or v == "":
                continue
            out[base] = v
        return out

    return wmap, getter


def form_values(wmap) -> dict:
    """从控件字典提取当前值 (False 视为关闭, 但保留在输出中)."""
    return {k: w.value for k, w in wmap.items()}

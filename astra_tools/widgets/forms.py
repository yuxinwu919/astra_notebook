"""元数据驱动的参数表单生成器 (ipywidgets).

从 astra_tools.deck.metadata 的参数元数据库 (类型/单位/默认值/描述)
自动渲染一个 namelist 的参数表单; 用户无需写代码即可填参数。
"""

from __future__ import annotations

import ipywidgets as widgets
import warnings as _warnings
from IPython.display import display

from ..deck.metadata import load


def _parse_fortran_number(s, integer: bool = False):
    """解析 Fortran 记法数值字面量.

    处理手册/历史 deck 中的常见记法:
      * 千分位分隔: '100 000', '500,000'
      * 指数 D/d:   '1.0D-3', '1.0d+3'
      * 浮点形式的整数: '1.0D+3' (integer=True 时取整)
    解析失败抛 ValueError —— 调用方决定降级策略, 不得静默回退 0
    (2026-08 审计 R2-1-1: sig_clock=0 会与 Cathode=T 组合产生简并
    束团, Max_step=0 使 ASTRA 立即终止)。
    """
    cleaned = str(s).strip().replace(",", "").replace(" ", "").replace("\t", "")
    if not cleaned:
        raise ValueError("empty number literal: %r" % (s,))
    cleaned = cleaned.replace("D", "E").replace("d", "e")
    if integer:
        try:
            return int(cleaned)
        except ValueError:
            return int(float(cleaned))
    return float(cleaned)


def _base_name(name: str) -> str:
    """参数基名: 去掉数组后缀 '()', '( )' 与二维数组手册记法 '(,)'.

    rstrip("() ") 无法处理 'D1(,)' (括号内含逗号), 二维数组参数
    (Q_mult_a(,)/Ap_GR(,) 等) 此前会得到 'D1(,' 的畸形基名。
    """
    n = name.rstrip()
    if n.endswith("(,)"):
        return n[:-3].rstrip()
    return n.rstrip("() ")


def _split_form_tokens(v: str) -> list:
    """逗号分隔的表单 token; 括号内 (复数/二维数组字面量) 逗号不切."""
    out = []
    cur = []
    depth = 0
    in_quote = None
    for ch in v:
        if in_quote:
            cur.append(ch)
            if ch == in_quote:
                in_quote = None
        elif ch in ("'", '"'):
            in_quote = ch
            cur.append(ch)
        elif ch == "(":
            depth += 1
            cur.append(ch)
        elif ch == ")":
            depth = max(0, depth - 1)
            cur.append(ch)
        elif ch == "," and depth == 0:
            out.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    out.append("".join(cur).strip())
    return [t for t in out if t]


def _form_complex(tok: str):
    """表单复数 token: '(1.0, 2.0)' -> 1+2j; 纯数 -> float."""
    t = tok.strip()
    if t.startswith("(") and t.endswith(")"):
        parts = [p.strip() for p in t[1:-1].split(",")]
        if len(parts) == 2:
            try:
                return complex(float(parts[0]), float(parts[1]))
            except ValueError:
                pass
    try:
        return float(t)
    except ValueError:
        raise ValueError(
            "参数需要逗号分隔的数值或 (re, im) 复数, 收到: %r" % tok)


def _widget_for_param(entry, value=None):
    """按参数元数据类型构造控件."""
    name = entry["name"]
    ptype = entry["type"]
    desc = entry["desc"]
    unit = entry["unit"] or ""
    tip = desc[:120] + ("  [" + unit + "]" if unit else "")

    if "Logical" in ptype:
        default = False
        if isinstance(value, bool):
            default = value
        elif isinstance(value, str):
            # 只有明确的真值拼写才勾选 (字符串 "FALSE" 不得翻成 True)
            default = value.strip().upper() in ("T", "TRUE")
        elif (isinstance(entry["default"], str)
              and entry["default"].strip().upper() in ("T", "TRUE")):
            default = True
        w = widgets.Checkbox(value=default, description=name, tooltip=tip)
    elif "Integer" in ptype and "array" not in ptype.lower():
        src = value if value is not None else entry.get("default")
        if src in (None, ""):
            w = widgets.IntText(value=0, description=name, tooltip=tip)
        else:
            try:
                w = widgets.IntText(
                    value=_parse_fortran_number(src, integer=True),
                    description=name, tooltip=tip)
            except (ValueError, TypeError):
                _warnings.warn(
                    "参数 %s: 默认值 %r 无法解析为整数, 保留原文显示 "
                    "(请手动修改, 不会静默回退 0)" % (name, src),
                    UserWarning, stacklevel=2)
                w = widgets.Text(
                    value=str(src), description=name, tooltip=tip,
                    layout=widgets.Layout(width="70%"))
    elif "Real" in ptype and "array" not in ptype.lower():
        src = value if value is not None else entry.get("default")
        if src in (None, ""):
            w = widgets.FloatText(value=0.0, description=name, tooltip=tip)
        else:
            try:
                w = widgets.FloatText(
                    value=_parse_fortran_number(src),
                    description=name, tooltip=tip)
            except (ValueError, TypeError):
                _warnings.warn(
                    "参数 %s: 默认值 %r 无法解析为实数, 保留原文显示 "
                    "(请手动修改, 不会静默回退 0)" % (name, src),
                    UserWarning, stacklevel=2)
                w = widgets.Text(
                    value=str(src), description=name, tooltip=tip,
                    layout=widgets.Layout(width="70%"))
    else:
        # Character / 数组 / 其他 -> 文本框 (数组用逗号分隔)
        default = ""
        if value is not None:
            if isinstance(value, (list, tuple)):
                default = ", ".join(str(v) for v in value)
            else:
                default = str(value)
        elif entry["default"] not in (None, ""):
            default = str(entry["default"])
            # 数值数组默认值可能含手册括号注释 (如 C_numb 的
            # "1 (3 in ASTRA Vers. 1)"), 清洗掉括号部分再放入控件
            if "array" in ptype.lower() and "Character" not in ptype:
                default = default.split("(")[0].strip()
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
    ptypes = {}
    touched = set()
    only_lower = None
    if only is not None:
        only_lower = [_base_name(o).lower() for o in only]
    matched_only = set()
    for entry in meta["params"]:
        base = _base_name(entry["name"])
        if only_lower is not None and base.lower() not in only_lower:
            continue
        matched_only.add(base.lower())
        val = None
        if values:
            for k, v in values.items():
                if _base_name(k).lower() == base.lower():
                    val = v
                    break
        wmap[base] = _widget_for_param(entry, val)
        ptypes[base] = entry["type"]
        # 记录用户改动 (changed_only 模式只写用户动过的参数)
        wmap[base].observe(lambda change, b=base: touched.add(b), "value")

    if only is not None:
        for o in only:
            o_base = _base_name(o).lower()
            if o_base not in matched_only:
                _warnings.warn(
                    "namelist %s: 参数 %r 不在元数据库中 (请检查拼写或补元数据)"
                    % (namelist, o), UserWarning, stacklevel=2)

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
            ptype = ptypes[base]
            # 文本框 -> 按元数据类型转值 (数组参数: 逗号分隔的数值
            # 列表; 否则写出带引号字符串会被 ASTRA 拒绝)
            if isinstance(v, str) and "array" in ptype.lower():
                toks = _split_form_tokens(v)
                if "Character" in ptype:
                    out[base] = toks
                elif "Complex" in ptype:
                    # 复数数组 (DIPOLE D1( , )): '(1.0, 2.0)' -> 1+2j
                    out[base] = [_form_complex(t) for t in toks]
                else:
                    try:
                        if "Integer" in ptype:
                            out[base] = [int(float(t)) for t in toks]
                        elif "Logical" in ptype:
                            out[base] = [t.upper() in ("T", "TRUE") for t in toks]
                        else:
                            out[base] = [float(t) for t in toks]
                    except ValueError:
                        raise ValueError(
                            "参数 %s 需要逗号分隔的数值, 收到: %r" % (base, v))
            else:
                out[base] = v
        return out

    return wmap, getter


def form_values(wmap) -> dict:
    """从控件字典提取当前值 (False 视为关闭, 但保留在输出中)."""
    return {k: w.value for k, w in wmap.items()}

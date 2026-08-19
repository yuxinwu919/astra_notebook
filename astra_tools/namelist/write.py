"""Write ASTRA/Generator namelist blocks from Python dictionaries.

Formatting rules (audited against the ASTRA Manual V3.2 examples):
  * booleans -> T / F
  * strings are quoted automatically unless they already contain quotes
    (ASTRA requires quoted strings, e.g. FNAME='bunch.ini')
  * ints and floats -> plain values (12 significant digits for floats)
  * complex -> Fortran literal (re, im)  (DIPOLE D1( , ) etc.)
  * lists/tuples/ndarrays -> comma-separated values; 2-D arrays are
    flattened column-major (Fortran order)
  * None and empty containers are skipped (optional parameters)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import numpy as np


def _format_complex(value) -> str:
    """Fortran 复数字面量 (re, im) — 手册 D1( , ) 写法. """
    return "(%s, %s)" % (format(float(value.real), ".12g"),
                         format(float(value.imag), ".12g"))


def _format_value(value: Any) -> str:
    if isinstance(value, (bool, np.bool_)):
        return "T" if value else "F"
    if isinstance(value, str):
        s = value.strip()
        if not s:
            return "''"  # 空串 -> 空字符字面量 (与跳过策略一致, 不抛)
        if s.startswith(("'", '"')) and s.endswith(("'", '"')):
            return s  # already quoted
        # Fortran 惯例: 内嵌撇号双写转义, 保证 parse 往返对称
        return "'" + s.replace("'", "''") + "'"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return format(float(value), ".12g")
    if isinstance(value, (complex, np.complexfloating)):
        return _format_complex(value)
    if isinstance(value, (list, tuple, np.ndarray)):
        # Fortran 是列主序: 二维数组按 order='F' 展平后整数组赋值
        # (2026-08 审计 R2-1-5: Q_mult_a(6,2)/D_Gap(2,n) 等)。
        # 不规则嵌套 (Module(, ) 字符二维数组) 用 object dtype 兜底。
        try:
            arr = np.asarray(value)
        except ValueError:
            arr = np.asarray(value, dtype=object)
        if arr.ndim >= 2:
            arr = arr.flatten(order="F")
        else:
            arr = arr.flatten()
        parts = []
        for v in arr:
            if isinstance(v, (bool, np.bool_)):
                parts.append("T" if v else "F")
            elif isinstance(v, str):
                parts.append(_format_value(v))
            elif isinstance(v, (int, np.integer)):
                parts.append(str(int(v)))
            elif isinstance(v, (complex, np.complexfloating)):
                parts.append(_format_complex(v))
            else:
                parts.append(format(float(v), ".12g"))
        return ", ".join(parts)
    raise TypeError("unsupported namelist value type: " + type(value).__name__)


def write_namelist(
    namelist_name: str,
    params: dict,
    filepath: Optional[Path] = None,
) -> Optional[str]:
    """Convert a parameter dict to an ASTRA-style namelist block.

    Args:
        namelist_name: block name, e.g. 'NEWRUN', 'INPUT'.
        params: dict of parameter name -> value.
        filepath: if given, append the block to this file (atomic-ish).

    Returns:
        The formatted text when filepath is None, else None.
    """
    lines = ["&" + namelist_name]
    for key, value in params.items():
        if value is None:
            continue
        if isinstance(value, (list, tuple, np.ndarray)) and len(value) == 0:
            continue
        if isinstance(value, (list, tuple, np.ndarray)) and len(value) == 1:
            # 单元素数组带 (1) 索引写出, parse 端恒存列表, 往返保持数组
            lines.append("  " + key + "(1)=" + _format_value(value) + ",")
            continue
        if isinstance(value, str) and value.strip() == "":
            continue
        lines.append("  " + key + "=" + _format_value(value) + ",")
    lines.append(" /")
    text = "\n".join(lines) + "\n"

    if filepath is None:
        return text
    filepath = Path(filepath)
    with open(filepath, "a") as f:
        f.write(text)
    return None


def write_input_deck(
    blocks: dict,
    filepath,
    header: str = "",
) -> Path:
    """Write a complete input deck: {block_name: params} -> file.

    Args:
        blocks: ordered dict of {namelist name: params dict}.
        filepath: target file (overwritten).
        header: optional leading comment text.
    """
    filepath = Path(filepath)
    with open(filepath, "w") as f:
        if header:
            for line in header.splitlines():
                f.write("! " + line + "\n")
        for name, params in blocks.items():
            if isinstance(params, (list, tuple)):
                # 重复 namelist 块 (手册第 3 章: INPUT 块重复 N_add 次)
                # 逐块写出, 保证 parse->write->parse 往返不变
                # (2026-08 审计 R2-1-3)。
                for p in params:
                    f.write(write_namelist(name, p))
                    f.write("\n")
            else:
                f.write(write_namelist(name, params))
                f.write("\n")
    return filepath

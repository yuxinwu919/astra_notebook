"""Write ASTRA/Generator namelist blocks from Python dictionaries.

Formatting rules (audited against the ASTRA Manual V3.2 examples):
  * booleans -> T / F
  * strings are quoted automatically unless they already contain quotes
    (ASTRA requires quoted strings, e.g. FNAME='bunch.ini')
  * ints and floats -> plain values (12 significant digits for floats)
  * lists/tuples/ndarrays -> comma-separated values
  * None and empty containers are skipped (optional parameters)
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

import numpy as np


def _format_value(value: Any) -> str:
    if isinstance(value, (bool, np.bool_)):
        return "T" if value else "F"
    if isinstance(value, str):
        s = value.strip()
        if not s:
            raise ValueError("empty string values are not valid namelist entries")
        if s.startswith(("'", '"')) and s.endswith(("'", '"')):
            return s  # already quoted
        # Fortran 惯例: 内嵌撇号双写转义, 保证 parse 往返对称
        return "'" + s.replace("'", "''") + "'"
    if isinstance(value, (int, np.integer)):
        return str(int(value))
    if isinstance(value, (float, np.floating)):
        return format(float(value), ".12g")
    if isinstance(value, (list, tuple, np.ndarray)):
        arr = np.asarray(value).flatten()
        parts = []
        for v in arr:
            if isinstance(v, (bool, np.bool_)):
                parts.append("T" if v else "F")
            elif isinstance(v, str):
                parts.append(_format_value(v))
            elif isinstance(v, (int, np.integer)):
                parts.append(str(int(v)))
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
        if isinstance(value, str) and value == "":
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
            f.write(write_namelist(name, params))
            f.write("\n")
    return filepath

"""Parse namelist blocks back into Python dicts (best-effort).

Useful for inspecting existing ASTRA input decks inside notebooks.
Values are parsed as int/float/bool/str; arrays become lists of floats
or ints. Unknown tokens stay strings.
"""

from __future__ import annotations

import math
import re
import warnings
from pathlib import Path


def _strip_comment(line: str) -> str:
    """在引号外找第一个 '!' 截断; 引号内的 '!' 保留."""
    out = []
    in_quote = None
    for ch in line:
        if in_quote:
            out.append(ch)
            if ch == in_quote:
                in_quote = None
        elif ch in ("'", '"'):
            in_quote = ch
            out.append(ch)
        elif ch == "!":
            break
        else:
            out.append(ch)
    return "".join(out).rstrip()


def _split_top_level(text: str, sep: str = ",") -> list:
    """按 sep 切分, 引号内的 sep 不切."""
    parts = []
    cur = []
    in_quote = None
    for ch in text:
        if in_quote:
            cur.append(ch)
            if ch == in_quote:
                in_quote = None
        elif ch in ("'", '"'):
            in_quote = ch
            cur.append(ch)
        elif ch == sep:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append("".join(cur))
    return parts


def _parse_token(tok: str):
    tok = tok.strip()
    if not tok:
        return None
    if tok[0] in ("'", '"') and tok[-1] == tok[0]:
        inner = tok[1:-1]
        # Fortran 双写转义: 'a''b' / "a""b"
        if tok[0] == "'":
            return inner.replace("''", "'")
        return inner.replace('""', '"')
    # 逻辑量: T/F 及 Fortran 点号形式 .T. / .F. (ASTRA 手册写法);
    # 点号形式必须在此识别, 否则会落成字符串 ".F" 并在写出时被
    # 加引号 (ADD='.F'), generator 会把引号内逻辑量判为非法并
    # 整体回退默认值 (Error reading input parameters)。
    if tok.upper() in ("T", ".T.", ".T", ".TRUE.", "TRUE"):
        return True
    if tok.upper() in ("F", ".F.", ".F", ".FALSE.", "FALSE"):
        return False
    try:
        if re.fullmatch(r"[+-]?\d+", tok):
            return int(tok)
        v = float(tok.replace("D", "E").replace("d", "e"))
        if math.isnan(v) or math.isinf(v):
            # 批 3: nan/inf 静默透传会污染重写后的 deck, 显式告警
            warnings.warn(
                "namelist token %r parses to %s (deck 将被写出 nan/inf)"
                % (tok, v), UserWarning, stacklevel=2)
        return v
    except ValueError:
        return tok


def parse_namelists(text_or_path) -> dict:
    """Parse an input deck into {block_name: {param: value}}.

    Array values (comma-separated) become lists. Repeated array
    elements with the same name (e.g. MaxE(1), MaxE(2)) are collected
    into one list.
    """
    if isinstance(text_or_path, Path) or (
        isinstance(text_or_path, str) and "\n" not in text_or_path
        and text_or_path.endswith((".in", ".dat"))
    ):
        text = Path(text_or_path).read_text()
    else:
        text = text_or_path

    blocks: dict = {}
    current = None
    last_base = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("!"):
            continue
        if line.startswith("&"):
            current = line[1:].strip()
            blocks[current] = {}
            last_base = None
            continue
        if line == "/":
            current = None
            last_base = None
            continue
        if current is None:
            continue
        content = _strip_comment(line).rstrip(",").strip()
        if "=" not in content:
            # 续行: 数组值换行延续 (手册: 逗号或换行皆可作分隔),
            # 无 '=' 的行追加到上一个参数 (此前被静默丢弃)
            if last_base is not None and content:
                tokens = [t for t in _split_top_level(content) if t.strip()]
                if tokens:
                    parsed = [_parse_token(t) for t in tokens]
                    prev = blocks[current].get(last_base)
                    if prev is None:
                        prev = []
                    elif not isinstance(prev, list):
                        prev = [prev]
                    prev.extend(parsed)
                    blocks[current][last_base] = prev
            continue

        def _store(key_raw, value_raw):
            nonlocal last_base
            key = key_raw.strip()
            m = re.search(r"\(\s*(\d+)\s*\)$", key)
            indexed = bool(m)
            # strip the (index) suffix from array elements
            base = re.sub(r"\(\s*\d+\s*\)$", "", key)
            last_base = base
            tokens = [t for t in _split_top_level(value_raw.strip()) if t.strip()]
            if not tokens:
                blocks[current][base] = None
                return
            parsed = [_parse_token(t) for t in tokens]
            if indexed:
                # 按索引放置 (乱序声明也正确); 单索引多值从 idx 起连续
                idx = int(m.group(1))
                arr = blocks[current].get(base)
                if not isinstance(arr, list):
                    arr = []
                while len(arr) < idx - 1 + len(parsed):
                    arr.append(None)
                for k, v in enumerate(parsed):
                    arr[idx - 1 + k] = v
                blocks[current][base] = arr
            else:
                if base in blocks[current]:
                    prev = blocks[current][base]
                    if not isinstance(prev, list):
                        prev = [prev]
                    prev.extend(parsed if len(parsed) > 1 else [parsed[0]])
                    blocks[current][base] = prev
                else:
                    blocks[current][base] = parsed[0] if len(parsed) == 1 else parsed

        # Multiple assignments may share one line (with or without
        # array parens), 引号感知:
        #   MaxE(1)=10, MaxE(2)=20,
        #   ZSTART=0.0, ZSTOP=1.5,
        #   Head='a, b', RUN=1,
        parts = _split_top_level(content)
        segments = []
        buf = []
        for seg in parts:
            if buf and re.match(r"^\s*[A-Za-z_]\w*\s*(?:\([^)]*\))?\s*=", seg):
                segments.append(",".join(buf).strip())
                buf = []
            buf.append(seg)
        if buf:
            segments.append(",".join(buf).strip())
        if len(segments) > 1:
            for seg in segments:
                seg = seg.rstrip(",").strip()
                if "=" in seg:
                    k, _, v = seg.partition("=")
                    _store(k, v)
        else:
            key, _, value = content.partition("=")
            _store(key, value)
    return blocks

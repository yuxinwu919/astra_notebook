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


def get_ci(d, key: str, default=None):
    """大小写不敏感读取 namelist 块参数 (ASTRA 键大小写不敏感).

    元数据/表单写出的键与手册一致 (LBfield), 旧 deck 可能用 LBField;
    parse_namelists 保留原始大小写, 因此查询必须不敏感
    (2026-08 审计 R2-1-2: 大小写不匹配导致 None -> 螺线管动量
    告警/修正被绕过)。
    """
    if not isinstance(d, dict):
        return default
    low = key.lower()
    for k, v in d.items():
        if k.lower() == low:
            return v
    return default


def iter_namelist_blocks(blocks: dict, name: str):
    """按声明顺序 yield 一个 namelist 的所有块.

    单块时 blocks[name] 是 dict -> yield 一次; 重复块时是 list[dict]
    (2026-08 审计 R2-1-3: 手册第 3 章要求 INPUT 块可重复 n 次)。
    """
    b = blocks.get(name)
    if b is None:
        return
    if isinstance(b, dict):
        yield b
    else:
        yield from b


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
    """按 sep 切分; 引号内与括号 ( ) 内的 sep 不切.

    括号感知支持复数/二维数组字面量 (1.0, 0.6) 与 D1(1,1) 类双下标键名
    (2026-08 审计 R2-1-5: DIPOLE D1=(1.0,2.0) 此前被切成 '(1.0' 与 '2.0)')。
    """
    parts = []
    cur = []
    in_quote = None
    depth = 0
    for ch in text:
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
        elif ch == sep and depth == 0:
            parts.append("".join(cur))
            cur = []
        else:
            cur.append(ch)
    parts.append("".join(cur))
    return parts


def _split_inline_terminator(line: str):
    """在引号/括号外找第一个 '/' (namelist 行内终止符, 手册示例
    'H_max=0.001, / ! end'); 返回 (terminator 之前的内容, 是否终止)。

    带行尾注释的裸终止行 '/ ! end of NEWRUN' 也由这里统一处理
    (2026-08 审计 R2-1-4: '/' 此前被当 token 塞进参数数组)。
    """
    depth = 0
    in_quote = None
    for i, ch in enumerate(line):
        if in_quote:
            if ch == in_quote:
                in_quote = None
            continue
        if ch in ("'", '"'):
            in_quote = ch
        elif ch == "(":
            depth += 1
        elif ch == ")":
            depth = max(0, depth - 1)
        elif ch == "/" and depth == 0:
            return line[:i], True
    return line, False


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
    # 复数 (手册 DIPOLE D1( , ): 'Two numbers, enclosed by brackets
    # and separated by a comma'): (1.0, 0.6) -> 1+0.6j
    if tok.startswith("(") and tok.endswith(")"):
        cparts = [p.strip() for p in tok[1:-1].split(",")]
        if len(cparts) == 2:
            try:
                re_ = float(cparts[0].replace("D", "E").replace("d", "e"))
                im_ = float(cparts[1].replace("D", "E").replace("d", "e"))
                return complex(re_, im_)
            except ValueError:
                pass
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


def _apply_line(cur_block: dict, content: str, last_base):
    """处理一行参数内容 (返回更新后的 last_base)."""
    if "=" not in content:
        # 续行: 数组值换行延续 (手册: 逗号或换行皆可作分隔),
        # 无 '=' 的行追加到上一个参数 (此前被静默丢弃)
        if last_base is not None and content:
            tokens = [t for t in _split_top_level(content) if t.strip()]
            if tokens:
                parsed = [_parse_token(t) for t in tokens]
                prev = cur_block.get(last_base)
                if prev is None:
                    prev = []
                elif not isinstance(prev, list):
                    prev = [prev]
                prev.extend(parsed)
                cur_block[last_base] = prev
        return last_base

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
            cur_block[base] = None
            return
        parsed = [_parse_token(t) for t in tokens]
        if indexed:
            # 按索引放置 (乱序声明也正确); 单索引多值从 idx 起连续
            idx = int(m.group(1))
            arr = cur_block.get(base)
            if not isinstance(arr, list):
                arr = []
            while len(arr) < idx - 1 + len(parsed):
                arr.append(None)
            for k, v in enumerate(parsed):
                arr[idx - 1 + k] = v
            cur_block[base] = arr
        else:
            if base in cur_block:
                prev = cur_block[base]
                if not isinstance(prev, list):
                    prev = [prev]
                prev.extend(parsed if len(parsed) > 1 else [parsed[0]])
                cur_block[base] = prev
            else:
                cur_block[base] = parsed[0] if len(parsed) == 1 else parsed

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
    return last_base


def parse_namelists(text_or_path) -> dict:
    """Parse an input deck into {block_name: {param: value}}.

    Array values (comma-separated) become lists. Repeated array
    elements with the same name (e.g. MaxE(1), MaxE(2)) are collected
    into one list. A namelist repeated in the deck (manual ch. 3:
    INPUT blocks repeat N_add times for multi-distribution generation)
    becomes a list of per-block dicts; single blocks stay dicts
    (2026-08 审计 R2-1-3)。
    """
    if isinstance(text_or_path, Path) or (
        isinstance(text_or_path, str) and "\n" not in text_or_path
        and text_or_path.endswith((".in", ".dat"))
    ):
        text = Path(text_or_path).read_text()
    else:
        text = text_or_path

    blocks: dict = {}
    cur_block: dict = None
    last_base = None
    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            continue
        stripped = _strip_comment(line)
        if not stripped:
            continue
        if stripped.startswith("&"):
            # 块名剥行尾注释: '&NEWRUN ! main block' -> 'NEWRUN'
            # (2026-08 审计 R2-1-4)
            name = stripped[1:].strip()
            new_block = {}
            prev = blocks.get(name)
            if prev is None:
                blocks[name] = new_block
            elif isinstance(prev, list):
                prev.append(new_block)
            else:
                blocks[name] = [prev, new_block]
            cur_block = new_block
            last_base = None
            continue
        content, terminated = _split_inline_terminator(stripped)
        if terminated:
            if cur_block is not None and content.strip():
                last_base = _apply_line(
                    cur_block, content.strip().rstrip(",").strip(), last_base)
            cur_block = None
            last_base = None
            continue
        if cur_block is None:
            continue
        content = content.rstrip(",").strip()
        if not content:
            continue
        last_base = _apply_line(cur_block, content, last_base)
    return blocks

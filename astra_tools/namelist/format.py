"""Namelist input-file formatter (normalizes ASTRA/Generator input decks).

Fixes the legacy format_input.py bugs:
  * '=' signs inside quoted strings are preserved
    (Head=' a = b ' stays intact)
  * trailing comments are recognized and no comma is appended after them
  * comments, block markers (&NAME, /) and '...' separators are preserved
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Optional


def _split_comment(line: str):
    """Split a line into (content, comment) outside of quotes."""
    in_quote = None
    for i, ch in enumerate(line):
        if in_quote:
            if ch == in_quote:
                in_quote = None
        elif ch in ("'", '"'):
            in_quote = ch
        elif ch == "!":
            return line[:i], line[i:]
    return line, ""


def _find_terminator(line: str):
    """在引号/括号外找第一个 '/' (namelist 行内终止符) 的位置.

    2026-08 审计 R2-1-4: 'H_max=0.001, / ! end' 的 '/' 不得被当作
    参数值保留 (否则回写 deck 时把 '/' 塞进参数数组)。
    """
    depth = 0
    in_quote = None
    for i, ch in enumerate(line):
        if in_quote:
            if ch == in_quote:
                in_quote = None
        elif ch in ("'", '"'):
            in_quote = ch
        elif ch == "(":
            depth += 1
        elif ch == ")":
            depth = max(0, depth - 1)
        elif ch == "/" and depth == 0:
            return i
    return None


def format_namelist_line(line: str) -> str:
    """Format one parameter line: 'key = value' -> 'key=value,'.

    Quotes are respected; a trailing comment is preserved without a
    comma being added after it.
    """
    stripped = line.strip()
    if not stripped:
        return ""
    if stripped.startswith("!"):
        return "  " + stripped
    if stripped.startswith("&"):
        return stripped
    if _split_comment(stripped)[0].strip() == "/":
        # 块终止行 (含尾注): '/ ! end' 保留注释
        # (2026-08 审计 R2-1-4: 此前注释被丢弃)
        _, term_comment = _split_comment(stripped)
        if term_comment.strip():
            return " / " + term_comment.strip()
        return " /"

    content, comment = _split_comment(stripped)
    content = content.strip()
    if not content:
        return "  " + comment if comment else ""

    # key = value -> key=value (only the first '=' outside quotes)
    if "=" in content:
        key, _, value = content.partition("=")
        key = key.rstrip()
        value = value.lstrip()
    else:
        key, value = content, ""
    if not value:
        return "  " + content

    formatted = key + "=" + value
    # remove spaces before '(' in array elements outside quotes:
    #   WK_Z (1) -> WK_Z(1)
    if "'" not in formatted and '"' not in formatted:
        formatted = re.sub(r"\s+\(", "(", formatted)

    if not formatted.endswith(","):
        formatted += ","
    if comment:
        formatted = formatted + " " + comment
    return "  " + formatted


def format_input_text(content: str) -> str:
    """Format a complete input deck.

    Normalizes line endings, indentation (2 spaces), trailing spaces,
    key=value spacing and commas; preserves comments, block markers and
    '...' separators.
    """
    content = content.replace("\r\n", "\n").replace("\r", "\n")
    result: list[str] = []
    in_namelist = False
    prev_blank = False

    for raw in content.split("\n"):
        line = raw.rstrip()
        stripped = line.strip()

        if stripped.startswith("&") and not stripped.startswith("&&"):
            in_namelist = True
            if result and result[-1] != "":
                result.append("")
            result.append(stripped)
            prev_blank = False
            continue
        if in_namelist:
            # 块终止: 裸 '/'、带尾注 '/ ! comment'、或行内
            # 'H_max=0.001, / ! end' (2026-08 审计 R2-1-4)。
            content, comment = _split_comment(stripped)
            term_idx = _find_terminator(content)
            if term_idx is not None:
                part = content[:term_idx].strip().rstrip(",").strip()
                if part:
                    formatted = format_namelist_line(part)
                    if formatted:
                        result.append(formatted)
                        prev_blank = False
                term = " /"
                if comment.strip():
                    term += " " + comment.strip()
                result.append(term)
                result.append("")
                in_namelist = False
                prev_blank = True
                continue
            formatted = format_namelist_line(stripped)
            if formatted:
                result.append(formatted)
                prev_blank = False
            continue
        # outside namelists
        if stripped:
            if stripped == "...":
                if result and result[-1] != "":
                    result.append("")
                result.append(stripped)
                result.append("")
                prev_blank = True
            else:
                result.append(stripped)
                prev_blank = False
        elif not prev_blank and result:
            result.append("")
            prev_blank = True

    while result and result[-1] == "":
        result.pop()
    return "\n".join(result) + "\n"


def format_input_file(
    filepath: Path,
    output: Optional[Path] = None,
    check_only: bool = False,
) -> bool:
    """Format a file in place (or to `output`).

    Returns True when the file was already normalized (no changes).
    """
    filepath = Path(filepath)
    try:
        original = filepath.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        return True  # not a text file
    formatted = format_input_text(original)
    if original == formatted:
        return True
    if check_only:
        return False
    target = Path(output) if output else filepath
    target.write_text(formatted, encoding="utf-8")
    return True

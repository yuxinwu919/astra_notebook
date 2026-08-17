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
    if stripped.startswith("&") or stripped == "/":
        return stripped

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
            # 块终止行可为裸 '/' 或带尾注的 '/ ! comment'
            if _split_comment(stripped)[0].strip() == "/":
                result.append(" /")
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

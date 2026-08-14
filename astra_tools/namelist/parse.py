"""Parse namelist blocks back into Python dicts (best-effort).

Useful for inspecting existing ASTRA input decks inside notebooks.
Values are parsed as int/float/bool/str; arrays become lists of floats
or ints. Unknown tokens stay strings.
"""

from __future__ import annotations

import re
from pathlib import Path


def _parse_token(tok: str):
    tok = tok.strip()
    if not tok:
        return None
    if tok[0] in ("'", '"') and tok[-1] == tok[0]:
        return tok[1:-1]
    if tok in ("T", ".TRUE.", "t", "true"):
        return True
    if tok in ("F", ".FALSE.", "f", "false"):
        return False
    try:
        if re.fullmatch(r"[+-]?\d+", tok):
            return int(tok)
        return float(tok.replace("D", "E").replace("d", "e"))
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
    for raw in text.splitlines():
        line = raw.strip()
        if not line or line.startswith("!"):
            continue
        if line.startswith("&"):
            current = line[1:].strip()
            blocks[current] = {}
            continue
        if line == "/":
            current = None
            continue
        if current is None:
            continue
        content = line.split("!", 1)[0].rstrip(",").strip()
        if "=" not in content:
            continue

        def _store(key_raw, value_raw):
            key = key_raw.strip()
            # strip the (index) suffix from array elements
            base = re.sub(r"\(\s*\d+\s*\)$", "", key)
            tokens = [t for t in re.split(r",\s*", value_raw.strip()) if t]
            if not tokens:
                blocks[current][base] = None
                return
            parsed = [_parse_token(t) for t in tokens]
            if base in blocks[current]:
                prev = blocks[current][base]
                if not isinstance(prev, list):
                    prev = [prev]
                prev.extend(parsed if len(parsed) > 1 else [parsed[0]])
                blocks[current][base] = prev
            else:
                blocks[current][base] = parsed[0] if len(parsed) == 1 else parsed

        # Multiple assignments may share one line:
        #   MaxE(1)=10, MaxE(2)=20,
        segments = re.split(r",\s*(?=[A-Za-z_]\w*\s*\()", content)
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

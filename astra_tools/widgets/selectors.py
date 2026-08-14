"""文件/run/z 位置三联动选择器 (ipywidgets).

ASTRA 输出命名: <stem>.<TYPE>.<run>, 相空间 <stem>.<NNNN>.<run>。
"""

from __future__ import annotations

from pathlib import Path

import ipywidgets as widgets

_TYPES = ["Xemit", "Yemit", "Zemit", "Cemit", "Sigma", "ref", "Log", "LandF"]


def discover_sim_runs(sim_dir):
    """扫描目录, 按 (stem, run) 分组输出文件.

    Returns:
        {stem: {run: {type: path, ...}, ...}, ...}
    """
    sim_dir = Path(sim_dir)
    found: dict = {}
    for f in sorted(sim_dir.iterdir()):
        if not f.is_file():
            continue
        parts = f.name.split(".")
        if len(parts) != 3 or not parts[2].isdigit():
            continue
        stem, typ, run = parts
        if typ in _TYPES or typ.isdigit():
            found.setdefault(stem, {}).setdefault(run, {})[typ] = f
    return found


def run_selector(sim_dir, stem: str = "", run: str = ""):
    """run 下拉选择器 (自动发现目录里的输出).

    Returns (selector, refresh_fn)。
    """
    runs = discover_sim_runs(sim_dir)
    stems = sorted(runs)
    if stem and stem not in stems:
        stems = [stem] + stems
    sel = widgets.Dropdown(
        options=stems,
        value=stem if stem in stems else (stems[0] if stems else None),
        description="run:",
        layout=widgets.Layout(width="auto"),
    )

    def refresh(new_dir=None):
        nonlocal runs
        if new_dir is not None:
            runs = discover_sim_runs(new_dir)
        st = sorted(runs)
        if not st:
            return
        sel.options = st
        if sel.value not in st:
            sel.value = st[0]

    return sel, refresh


def phase_selector(phase_files):
    """相空间 z 位置选择器.

    Args:
        phase_files: 相空间文件 Path 列表 (如 astra.0100.001)。
    """
    options = []
    for f in sorted(phase_files):
        tag = f.name.split(".")[1]
        try:
            z_cm = int(tag)
            label = "z = %.3f m  (%s)" % (z_cm / 100.0, f.name)
        except ValueError:
            label = f.name
        options.append((label, f))
    sel = widgets.Select(
        options=options,
        value=options[-1][1] if options else None,
        description="z 位置:",
        rows=min(8, max(len(options), 1)),
        layout=widgets.Layout(width="90%"),
    )
    return sel

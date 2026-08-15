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
        if typ in _TYPES or typ.lstrip("-").isdigit():
            found.setdefault(stem, {}).setdefault(run, {})[typ] = f
    return found


def run_selector(sim_dir, stem: str = "", run: str = ""):
    """算例(stem) 下拉选择器 (自动发现目录里的输出).

    选项是输出文件的项目名 stem (如 Example、Wake); 同 stem 的
    不同 run 由 phase/z 选择器再细分。run 参数为兼容保留。

    Returns (selector, refresh_fn)。
    """
    runs = discover_sim_runs(sim_dir)
    stems = sorted(runs)
    if stem and stem not in stems:
        stems = [stem] + stems
    sel = widgets.Dropdown(
        options=stems,
        value=stem if stem in stems else (stems[0] if stems else None),
        description="算例:",
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


def _phase_label(f):
    """相空间文件 -> 显示标签 (z 位置)。"""
    try:
        z_cm = int(f.name.split(".")[1])
        return "z = %.3f m  (%s)" % (z_cm / 100.0, f.name)
    except ValueError:
        return f.name


class PhaseStepper(widgets.VBox):
    """z 位置步进器 (对应 postpro 的步进功能)。

    滑块 (0..N-1) 是状态源; ⏮◀▶⏭ 按钮切换; 底部标签显示当前 z。
    属性: .index (IntSlider), .path (当前相空间文件 Path)。配合
    ipywidgets.interactive_output 可步进时自动刷新统计/相空间
    (见 03_postpro.ipynb)。
    """

    def __init__(self, phase_files, start: int = -1):
        super().__init__()
        self.files = sorted(phase_files)
        self.labels = [_phase_label(f) for f in self.files]
        n = max(len(self.files), 1)
        i0 = max(min(start if start >= 0 else n - 1, n - 1), 0)
        self.index = widgets.IntSlider(
            min=0, max=n - 1, value=i0, description="", readout=False,
            continuous_update=False,
            layout=widgets.Layout(width="220px"),
        )
        self.lab = widgets.Label(
            value=self.labels[i0] if self.files else "无相空间文件")
        btn_first = widgets.Button(description="⏮", tooltip="第一个 z 位置",
                                   layout=widgets.Layout(width="40px"))
        btn_prev = widgets.Button(description="◀", tooltip="上一个 z 位置 (步进)",
                                  layout=widgets.Layout(width="40px"))
        btn_next = widgets.Button(description="▶", tooltip="下一个 z 位置 (步进)",
                                  layout=widgets.Layout(width="40px"))
        btn_last = widgets.Button(description="⏭", tooltip="最后一个 z 位置",
                                  layout=widgets.Layout(width="40px"))

        def _go(delta):
            if self.files:
                self.index.value = min(max(self.index.value + delta, 0), n - 1)

        btn_first.on_click(lambda _b: setattr(self.index, "value", 0))
        btn_prev.on_click(lambda _b: _go(-1))
        btn_next.on_click(lambda _b: _go(+1))
        btn_last.on_click(lambda _b: setattr(self.index, "value", n - 1))

        def _sync(_change):
            if self.files:
                self.lab.value = self.labels[self.index.value]

        self.index.observe(_sync, "value")
        self.children = [
            widgets.HBox([btn_first, btn_prev, self.index, btn_next, btn_last]),
            self.lab,
        ]

    @property
    def label(self) -> str:
        """当前显示标签 (z 位置)。"""
        return self.lab.value

    @property
    def path(self):
        """当前相空间文件路径 (无文件时为 None)。"""
        if not self.files:
            return None
        return self.files[self.index.value]


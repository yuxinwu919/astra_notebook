"""文件/run/z 位置三联动选择器 (ipywidgets).

ASTRA 输出命名: <stem>.<TYPE>.<run>, 相空间 <stem>.<NNNN>.<run>。
"""

from __future__ import annotations

from pathlib import Path

import ipywidgets as widgets
import numpy as np

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
    options = [(_phase_label(f), f) for f in sorted(phase_files)]

    sel = widgets.Select(
        options=options,
        value=options[-1][1] if options else None,
        description="z 位置:",
        rows=min(8, max(len(options), 1)),
        layout=widgets.Layout(width="90%"),
    )
    return sel


def _phase_label(f):
    """相空间文件 -> 显示标签 (z 位置)。

    批 6: 优先读文件首行绝对 z (ASCII 首行=参考粒子, 第 3 列);
    二进制回退文件名启发 (cm: 0150=1.5m 与 mm: 1500=1.5m 双约定,
    此前一律按 cm 除 100, mm 命名的 dump 会显示 10 倍大的 z)。
    """
    try:
        with open(f, "rb") as fh:
            head = fh.read(4096)
        if b"\x00" not in head:
            rows = [ln for ln in head.decode("ascii", "ignore").splitlines()
                    if ln.strip() and not ln.strip().startswith(("#", "!"))]
            if rows:
                z_m = float(rows[0].split()[2])
                return "z = %.4f m  (%s)" % (z_m, f.name)
    except Exception:
        pass
    try:
        v = int(f.name.split(".")[1])
        z_m = v / 1000.0 if v >= 1000 else v / 100.0
        return "z = %.4f m  (%s)" % (z_m, f.name)
    except (ValueError, IndexError):
        return f.name


class PhaseStepper(widgets.VBox):
    """z 位置步进器 (对应 postpro 的步进功能)。

    滑块 (0..N-1) 是状态源; ◀◀ ◀ ▶ ▶▶ 按钮切换; 底部标签显示当前 z。
    属性: .index (IntSlider), .path (当前相空间文件 Path)。配合
    ipywidgets.interactive_output 可步进时自动刷新统计/相空间
    (见 postpro.ipynb)。
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
        btn_first = widgets.Button(description="◀◀", tooltip="第一个 z 位置",
                                   layout=widgets.Layout(width="40px"))
        btn_prev = widgets.Button(description="◀", tooltip="上一个 z 位置 (步进)",
                                  layout=widgets.Layout(width="40px"))
        btn_next = widgets.Button(description="▶", tooltip="下一个 z 位置 (步进)",
                                  layout=widgets.Layout(width="40px"))
        btn_last = widgets.Button(description="▶▶", tooltip="最后一个 z 位置",
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


_PARAM_NAMES = ["x", "y", "z", "px", "py", "pz", "clock", "t",
                "xp", "yp", "dp/p", "E_kin"]


class PhaseSpaceParamSelector(widgets.VBox):
    """任意参数对相空间控件 (postpro 5.6.2 菜单 2).

    两个参数下拉 + 加投影/减线性相关/状态着色 复选框; 绘图按钮;
    叠加管理 (保存层 / 叠加图 / 清空)。配合 Output 使用。
    """

    def __init__(self, dist, bz_on_axis_T: float = 0.0):
        super().__init__()
        from ..plot.arbitrary_phase_space import OverlayManager
        self.dist = dist
        self.bz = bz_on_axis_T
        self.overlay = OverlayManager()
        self.x_sel = widgets.Dropdown(options=_PARAM_NAMES, value="x",
                                      description="X", layout=widgets.Layout(width="130px"))
        self.y_sel = widgets.Dropdown(options=_PARAM_NAMES, value="xp",
                                      description="Y", layout=widgets.Layout(width="130px"))
        self.sub_corr = widgets.Checkbox(value=False, description="减线性相关")
        self.add_proj = widgets.Checkbox(value=False, description="加投影")
        self.color_status = widgets.Checkbox(value=False, description="状态着色")
        self.btn_plot = widgets.Button(description="绘图", layout=widgets.Layout(width="70px"))
        self.btn_add = widgets.Button(description="加入叠加", layout=widgets.Layout(width="90px"))
        self.btn_overlay = widgets.Button(description="叠加图", layout=widgets.Layout(width="70px"))
        self.btn_clear = widgets.Button(description="清空", layout=widgets.Layout(width="60px"))
        self.lab = widgets.Label(value="叠加层数: 0")

        def _plot(_b=None):
            from ..plot.arbitrary_phase_space import plot_arbitrary
            return plot_arbitrary(
                self.dist, self.x_sel.value, self.y_sel.value,
                subtract_corr=self.sub_corr.value, add_proj=self.add_proj.value,
                color_by_status=self.color_status.value, bz_on_axis_T=self.bz)

        def _add(_b=None):
            self.overlay.add(self.dist, self.x_sel.value, self.y_sel.value,
                             bz_on_axis_T=self.bz)
            self.lab.value = "叠加层数: %d" % self.overlay.count

        def _show(_b=None):
            return self.overlay.plot()

        self.btn_plot.on_click(_plot)
        self.btn_add.on_click(_add)
        self.btn_overlay.on_click(_show)
        self.btn_clear.on_click(lambda _b: (self.overlay.clear(),
                                            setattr(self.lab, "value", "叠加层数: 0")))
        self.children = [
            widgets.HBox([self.x_sel, self.y_sel]),
            widgets.HBox([self.sub_corr, self.add_proj, self.color_status]),
            widgets.HBox([self.btn_plot, self.btn_add, self.btn_overlay,
                          self.btn_clear, self.lab]),
        ]


class CutControls(widgets.VBox):
    """相空间切割控件 (postpro 5.6.4, 滑块替代鼠标).

    x/y/z/E 四组窗口滑块 + 应用/撤销。保留原分布副本用于撤销。
    属性: .dist (当前分布), .original。
    """

    def __init__(self, dist):
        super().__init__()
        from IPython.display import display
        self.original = dist
        self.dist = dist
        self._applied = False

        def _mk(param, label, lo, hi, unit, n=200):
            return widgets.FloatRangeSlider(
                min=lo, max=hi, value=(lo, hi), step=(hi - lo) / n,
                description=label, continuous_update=False,
                layout=widgets.Layout(width="420px"))

        d = dist.filter_active()
        self.sliders = {
            "x": _mk("x", "x [mm]", float(d.x.min() * 1e3), float(d.x.max() * 1e3), "mm"),
            "y": _mk("y", "y [mm]", float(d.y.min() * 1e3), float(d.y.max() * 1e3), "mm"),
            "z": _mk("z", "z [mm]", float(d.z.min() * 1e3), float(d.z.max() * 1e3), "mm"),
            "E": _mk("E", "E [MeV]", float(
                _e(d).min() * 1e-6), float(_e(d).max() * 1e-6), "MeV"),
        }
        self.btn_apply = widgets.Button(description="应用切割")
        self.btn_reject = widgets.Button(description="撤销")
        self.lab = widgets.Label(value="")

        def _apply(_b=None):
            from ..analysis.cuts import cut_distribution
            from ..constants import MEV_TO_EV
            lo, hi = self.sliders["x"].value
            xr = (lo * 1e-3, hi * 1e-3)
            lo, hi = self.sliders["y"].value
            yr = (lo * 1e-3, hi * 1e-3)
            lo, hi = self.sliders["z"].value
            zr = (lo * 1e-3, hi * 1e-3)
            lo, hi = self.sliders["E"].value
            er = (lo * MEV_TO_EV, hi * MEV_TO_EV)
            self.dist, _mask = cut_distribution(
                self.original, x_range=xr, y_range=yr, z_range=zr, e_range=er)
            self._applied = True
            self.lab.value = "已应用切割 (active: %d/%d)" % (
                self.dist.n_active, self.original.n_active)

        def _reject(_b=None):
            if self._applied:
                self.dist = self.original
                self._applied = False
                self.lab.value = "已撤销, 恢复原分布"

        self.btn_apply.on_click(_apply)
        self.btn_reject.on_click(_reject)
        self.children = [
            self.sliders["x"], self.sliders["y"],
            self.sliders["z"], self.sliders["E"],
            widgets.HBox([self.btn_apply, self.btn_reject, self.lab]),
        ]


def _e(dist):
    """active 粒子动能 [eV] (CutControls 用)."""
    from ..constants import kinetic_energy_from_momentum
    m = dist.active
    return np.asarray(kinetic_energy_from_momentum(dist.pz[m]), dtype=float)


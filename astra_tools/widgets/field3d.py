"""3D 场图 2D 交互剖面 (原 fieldplot 菜单 2 的现代版: 选平面 + 滑块选位置).

需要 ipywidgets; 无 Jupyter 环境下可导入, 组件仅在 notebook 中显示。
"""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt
from IPython.display import display

from ..io.field_map import read_3d_field_map_components
from ..plot.field_plots import plot_3d_field_slices, resolve_3d_plane, plane_fixed_axis


def interact_3d_field_slices(base, kind="heatmap", plane="auto",
                             auto_render=True):
    """交互式 3D 场图剖面: 平面单选 + 位置滑块步进 + 画法切换.

    Args:
        base: 场图主名 (3D_Dipole) 或任一分量文件 (3D_Dipole.by)。
        kind: 初始画法 'heatmap' (|F| 色块) | 'contour' (填充等值线)。
        plane: 初始平面 'auto'/'xy'/'xz'/'yz'。
        auto_render: 构造时是否立即渲染 (测试传 False)。

    Returns:
        ipywidgets.VBox — notebook 中 display 后即可交互。
    """
    import ipywidgets as widgets

    field = read_3d_field_map_components(base)
    start_plane = resolve_3d_plane(field, plane, mode="plane")

    def fixed_grid(pl):
        fixed_axis = plane_fixed_axis(pl)
        return fixed_axis, getattr(field, fixed_axis)

    w_kind = widgets.ToggleButtons(
        options=["heatmap", "contour"], value=kind, description="画法:")
    w_plane = widgets.ToggleButtons(
        options=["xy", "xz", "yz"], value=start_plane, description="平面:")

    fixed_axis, fixed = fixed_grid(start_plane)
    w_slider = widgets.FloatSlider(
        value=(float(fixed.min()) + float(fixed.max())) * 0.5 * 1e3,
        min=float(fixed.min() * 1e3),
        max=float(fixed.max() * 1e3),
        step=float(np.median(np.diff(fixed))) * 1e3,
        description="%s [mm]:" % fixed_axis,
        layout=widgets.Layout(width="85%"), continuous_update=False)
    out = widgets.Output()
    updating_plane = False

    def render():
        pl = w_plane.value
        fig = plot_3d_field_slices(field, plane=pl, kind=w_kind.value,
                                   position=float(w_slider.value) * 1e-3)
        with out:
            out.clear_output(wait=True)
            display(fig)
        # 立即关闭 figure: 否则 matplotlib inline 的 flush_figures (post_execute)
        # 会在 cell 执行结束时把仍打开的图再次显示成新的 cell 输出 (重复图片)。
        plt.close(fig)

    def on_plane(_):
        nonlocal updating_plane
        fa, fg = fixed_grid(w_plane.value)
        updating_plane = True
        try:
            w_slider.min = float(fg.min() * 1e3)
            w_slider.max = float(fg.max() * 1e3)
            w_slider.step = float(np.median(np.diff(fg))) * 1e3
            w_slider.description = "%s [mm]:" % fa
            w_slider.value = (w_slider.min + w_slider.max) * 0.5
        finally:
            updating_plane = False
        render()

    def on_slider(_):
        if not updating_plane:
            render()

    w_plane.observe(on_plane, names="value")
    w_slider.observe(on_slider, names="value")
    w_kind.observe(lambda *_: render(), names="value")
    if auto_render:
        render()
    return widgets.VBox([w_kind, w_plane, w_slider, out])


def interact_3d_field_map(base, plane="auto", kind="heatmap",
                          auto_render=True):
    """交互式场图 2D 剖面: 平面切换 + 位置滑块步进 + 画法切换。

    拖动滑块只刷新当前 2D 切片，不再绘制 3D 全景或面内箭头。

    Args:
        base: 场图主名 (3D_Dipole) 或任一分量文件 (3D_Dipole.by)。
        plane: 初始平面 'auto'/'xy'/'xz'/'yz'。
        kind: 2D 剖面初始画法 'heatmap' (|F| 色块) | 'contour'。
        auto_render: 构造时是否立即渲染 (测试传 False)。

    Returns:
        ipywidgets.VBox — notebook 中 display 后即可交互。
    """
    import ipywidgets as widgets

    field = read_3d_field_map_components(base)
    start_plane = resolve_3d_plane(field, plane, mode="plane")

    def fixed_grid(pl):
        fixed_axis = plane_fixed_axis(pl)
        return fixed_axis, getattr(field, fixed_axis)

    w_kind = widgets.ToggleButtons(options=["heatmap", "contour"],
                                   value=kind, description="画法:")
    w_plane = widgets.ToggleButtons(options=["xy", "xz", "yz"],
                                    value=start_plane, description="平面:")
    fixed_axis, fixed = fixed_grid(start_plane)
    w_slider = widgets.FloatSlider(
        value=(float(fixed.min()) + float(fixed.max())) * 0.5 * 1e3,
        min=float(fixed.min() * 1e3),
        max=float(fixed.max() * 1e3),
        step=float(np.median(np.diff(fixed))) * 1e3,
        description="%s [mm]:" % fixed_axis,
        layout=widgets.Layout(width="85%"),
        continuous_update=False)
    out = widgets.Output()
    updating_plane = False

    def render():
        pl = w_plane.value
        pos_m = float(w_slider.value) * 1e-3
        fig2d = plot_3d_field_slices(field, plane=pl, kind=w_kind.value,
                                     position=pos_m, figsize=(7.5, 3.6))
        with out:
            out.clear_output(wait=True)
            display(fig2d)
        # 立即关闭 figure: 否则 matplotlib inline 的 flush_figures (post_execute)
        # 会在 cell 执行结束时把仍打开的图再次显示成新的 cell 输出 (重复图片)。
        plt.close(fig2d)

    def on_plane(_):
        nonlocal updating_plane
        fa, fg = fixed_grid(w_plane.value)
        updating_plane = True
        try:
            w_slider.min = float(fg.min() * 1e3)
            w_slider.max = float(fg.max() * 1e3)
            w_slider.step = float(np.median(np.diff(fg))) * 1e3
            w_slider.description = "%s [mm]:" % fa
            w_slider.value = (w_slider.min + w_slider.max) * 0.5
        finally:
            updating_plane = False
        render()

    def on_slider(_):
        if not updating_plane:
            render()

    w_plane.observe(on_plane, names="value")
    w_slider.observe(on_slider, names="value")
    w_kind.observe(lambda *_: render(), names="value")
    if auto_render:
        render()
    return widgets.VBox([w_kind, w_plane, w_slider, out])

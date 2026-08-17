"""3D 场图交互剖面 (原 fieldplot 菜单 2 的现代版: 选平面 + 滑块选位置).

需要 ipywidgets; 无 Jupyter 环境下可导入, 组件仅在 notebook 中显示。
"""

from __future__ import annotations

import numpy as np
from IPython.display import display

from ..io.field_map import read_3d_field_map_components
from ..plot.field_plots import (plot_3d_field_slices, plot_3d_field_contour,
                                _coerce_plane, _PLANE_SPEC)


def interact_3d_field_slices(base, view="vector_slices", plane="auto",
                             auto_render=True):
    """交互式 3D 场图剖面: 平面单选 + 位置滑块 + 视图切换.

    Args:
        base: 场图主名 (3D_Dipole) 或任一分量文件 (3D_Dipole.by)。
        view: 初始视图 'vector_slices' | 'contour'。
        plane: 初始平面 'auto'/'xy'/'xz'/'yz'。
        auto_render: 构造时是否立即渲染 (测试传 False)。

    Returns:
        ipywidgets.VBox — notebook 中 display 后即可交互。
    """
    import ipywidgets as widgets

    field = read_3d_field_map_components(base)
    start_plane = _coerce_plane(field, plane, None, mode="plane")

    def fixed_grid(pl):
        fixed_axis = "xyz"[_PLANE_SPEC[pl][0]]
        return fixed_axis, getattr(field, fixed_axis)

    w_view = widgets.ToggleButtons(
        options=["vector_slices", "contour"], value=view, description="视图:")
    w_plane = widgets.ToggleButtons(
        options=["xy", "xz", "yz"], value=start_plane, description="平面:")

    fixed_axis, fixed = fixed_grid(start_plane)
    w_slider = widgets.FloatSlider(
        value=(float(fixed.min()) + float(fixed.max())) * 0.5,
        min=float(fixed.min() * 1e3),
        max=float(fixed.max() * 1e3),
        step=float(np.median(np.diff(fixed))) * 1e3,
        description="%s [mm]:" % fixed_axis,
        layout=widgets.Layout(width="85%"))
    out = widgets.Output()

    def render():
        pl = w_plane.value
        with out:
            out.clear_output(wait=True)
            pos_m = float(w_slider.value) * 1e-3
            if w_view.value == "contour":
                fig = plot_3d_field_contour(field, plane=pl, position=pos_m)
            else:
                fig = plot_3d_field_slices(field, plane=pl, position=pos_m)
            display(fig)

    def on_plane(_):
        fa, fg = fixed_grid(w_plane.value)
        w_slider.min = float(fg.min() * 1e3)
        w_slider.max = float(fg.max() * 1e3)
        w_slider.step = float(np.median(np.diff(fg))) * 1e3
        w_slider.description = "%s [mm]:" % fa
        w_slider.value = (w_slider.min + w_slider.max) * 0.5
        render()

    w_plane.observe(on_plane, names="value")
    w_slider.observe(lambda *_: render(), names="value")
    w_view.observe(lambda *_: render(), names="value")
    if auto_render:
        render()
    return widgets.VBox([w_view, w_plane, w_slider, out])

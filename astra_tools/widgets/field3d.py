"""3D 场图交互剖面 (原 fieldplot 菜单 2 的现代版: 选平面 + 滑块选位置).

需要 ipywidgets; 无 Jupyter 环境下可导入, 组件仅在 notebook 中显示。
"""

from __future__ import annotations

import numpy as np
from IPython.display import display

from ..io.field_map import read_3d_field_map_components
from ..plot.field_plots import (plot_3d_field_slices, plot_3d_field_quiver,
                                plot_3d_field_stack, resolve_3d_plane,
                                plane_fixed_axis)


def interact_3d_field_slices(base, kind="heatmap", plane="auto",
                             auto_render=True):
    """交互式 3D 场图剖面: 平面单选 + 位置滑块步进 + 画法切换.

    Args:
        base: 场图主名 (3D_Dipole) 或任一分量文件 (3D_Dipole.by)。
        kind: 初始画法 'heatmap' (|F| 色块 + 箭头) | 'contour' (填充等值线)。
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
            fig = plot_3d_field_slices(field, plane=pl, kind=w_kind.value,
                                       position=pos_m)
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
    w_kind.observe(lambda *_: render(), names="value")
    if auto_render:
        render()
    return widgets.VBox([w_kind, w_plane, w_slider, out])


def interact_3d_field_quiver(base, plane="xy", auto_render=True):
    """交互式 3D 场图 2D quiver: 平面切换 + 位置滑块步进.

    选定显示平面 xy/xz/yz, 显示面内 2D quiver 剖面 (如 plane='xy' ->
    箭头 (U,V) = (fx, fy), 即 (Ex,Ey)/(Bx,By)); 滑块控制固定轴切片位置,
    拖动自动重绘。纯 2D quiver (无 |F| 底色)。

    Args:
        base: 场图主名 (3D_Dipole) 或任一分量文件 (3D_Dipole.by)。
        plane: 初始平面 'xy'/'xz'/'yz'。
        auto_render: 构造时是否立即渲染 (测试传 False)。

    Returns:
        ipywidgets.VBox — notebook 中 display 后即可交互。
    """
    import ipywidgets as widgets

    field = read_3d_field_map_components(base)
    if plane not in ("xy", "xz", "yz"):
        raise ValueError("plane 必须为 'xy'/'xz'/'yz': %r" % (plane,))

    w_plane = widgets.ToggleButtons(
        options=["xy", "xz", "yz"], value=plane, description="平面:")
    fixed_axis = plane_fixed_axis(plane)
    fixed = getattr(field, fixed_axis)
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
            fig = plot_3d_field_quiver(field, plane=pl, position=pos_m)
            display(fig)

    def on_plane(_):
        fa = plane_fixed_axis(w_plane.value)
        fg = getattr(field, fa)
        w_slider.min = float(fg.min() * 1e3)
        w_slider.max = float(fg.max() * 1e3)
        w_slider.step = float(np.median(np.diff(fg))) * 1e3
        w_slider.description = "%s [mm]:" % fa
        w_slider.value = (w_slider.min + w_slider.max) * 0.5
        render()

    w_plane.observe(on_plane, names="value")
    w_slider.observe(lambda *_: render(), names="value")
    if auto_render:
        render()
    return widgets.VBox([w_plane, w_slider, out])


def interact_3d_field_map(base, plane="auto", kind="heatmap", n_slices=7,
                          auto_render=True):
    """统一交互式 3D 场图: 平面切换 + 位置滑块步进 + 画法切换.

    上图为 3D 上下文 (半透明切片栈 + 当前切片高亮平面), 下图为 2D 剖面
    (当前切片的 heatmap/contour + 面内箭头)。拖动滑块同时刷新上下两图,
    把「3D 全景」和「切片步进」合成到一个组件里。

    Args:
        base: 场图主名 (3D_Dipole) 或任一分量文件 (3D_Dipole.by)。
        plane: 初始平面 'auto'/'xy'/'xz'/'yz'。
        kind: 2D 剖面初始画法 'heatmap' (|F| 色块+箭头) | 'contour'。
        n_slices: 3D 上下文里的参考切片层数。
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
        value=(float(fixed.min()) + float(fixed.max())) * 0.5,
        min=float(fixed.min() * 1e3),
        max=float(fixed.max() * 1e3),
        step=float(np.median(np.diff(fixed))) * 1e3,
        description="%s [mm]:" % fixed_axis,
        layout=widgets.Layout(width="85%"),
        continuous_update=False)
    out = widgets.Output()

    def render():
        pl = w_plane.value
        pos_m = float(w_slider.value) * 1e-3
        with out:
            out.clear_output(wait=True)
            fig3d = plot_3d_field_stack(field, plane=pl, n_slices=n_slices,
                                        highlight=pos_m, figsize=(7.5, 4.8))
            display(fig3d)
            fig2d = plot_3d_field_slices(field, plane=pl, kind=w_kind.value,
                                         position=pos_m, figsize=(7.5, 3.6))
            display(fig2d)

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
    w_kind.observe(lambda *_: render(), names="value")
    if auto_render:
        render()
    return widgets.VBox([w_kind, w_plane, w_slider, out])

"""postpro P2 遗留项测试 (slice 投影椭圆 2D, CP_ind 粒子索引着色).

手册 5.6.3 项 4 (投影 rms slice 椭圆) / 项 7 (投影切换) / 项 11 (减线性相关)
与 5.6 (Plot_steering.par CP_ind, Plot_mode=1)。
"""

from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]

from astra_tools.io import read_distribution, cp_index_colors
from astra_tools.plot.advanced_plots import plot_slice_ellipses_2d
from astra_tools.plot.phase_space import plot_phase_space

DATA = PROJECT_ROOT / "examples/Manual_Example"


@pytest.fixture(scope="module")
def dist():
    return read_distribution(DATA / "Example.0150.001")


def test_cp_index_colors_mapping():
    idx = np.array([1, 1, 2, 3, 1])
    colors = cp_index_colors(idx, {1: (1, 0, 0), 2: (0, 1, 0)})
    assert colors[0] == (1, 0, 0)
    assert colors[1] == (1, 0, 0)
    assert colors[2] == (0, 1, 0)
    assert colors[3] == "#0077BB"   # 未指定索引 -> fallback
    # 黑色 (0,0,0) = 不绘制 -> None
    colors2 = cp_index_colors(idx, {1: (0, 0, 0), 2: (1, 0, 0)})
    assert colors2[0] is None
    assert colors2[2] == (1, 0, 0)


def test_slice_ellipses_2d_renders(dist):
    import matplotlib.pyplot as plt
    fig = plot_slice_ellipses_2d(dist, n_slices=6, plane="xxp")
    try:
        assert "x [mm]" in fig.axes[0].get_xlabel()
        assert "mrad" in fig.axes[0].get_ylabel()
        # 散点 collection + n_slices 椭圆 line
        n_ell = sum(1 for ln in fig.axes[0].lines
                    if "slice" not in ln.get_label())
        assert n_ell >= 1
    finally:
        plt.close("all")


def test_slice_ellipses_2d_planes(dist):
    import matplotlib.pyplot as plt
    for plane in ("xxp", "yyp", "xyp", "yxp"):
        fig = plot_slice_ellipses_2d(dist, n_slices=5, plane=plane)
        try:
            assert any(p in fig.axes[0].get_ylabel()
                       for p in ("mrad", "mm"))
        finally:
            plt.close("all")


def test_slice_ellipses_2d_subtract_corr(dist):
    import matplotlib.pyplot as plt
    fig = plot_slice_ellipses_2d(dist, n_slices=6, subtract_corr=True)
    try:
        assert "corr removed" in fig.axes[0].get_title()
    finally:
        plt.close("all")


def test_phase_space_colors_array(dist):
    """逐粒子颜色数组传给 plot_phase_space 不抛异常."""
    import matplotlib.pyplot as plt
    colors = cp_index_colors(dist.index, {1: (1, 0, 0), 2: (0, 1, 0)})
    fig = plot_phase_space(dist, plane="x", colors=colors,
                           title="colored by CP_ind")
    try:
        assert len(fig.axes[0].collections) >= 1
    finally:
        plt.close("all")


def test_phase_space_colors_black_none(dist):
    """含黑色 (None) 的颜色数组: 调用方替换后仍可渲染."""
    import matplotlib.pyplot as plt
    colors = cp_index_colors(dist.index, {1: (0, 0, 0)})
    colors = np.array([c if c is not None else "#0077BB" for c in colors],
                      dtype=object)
    fig = plot_phase_space(dist, plane="x", colors=colors)
    try:
        assert len(fig.axes[0].collections) >= 1
    finally:
        plt.close("all")

"""第二轮审计 UI 层修复测试 (R2-2: notebook 前端数据流).

覆盖:
  * R2-2-1  phase_file_z_m: 相空间文件名 -> 束团 z [m] (cm/mm 双约定,
            复用 _phase_label 的解析逻辑); 非相位文件名 -> None。
            P3: _phase_label 接受 str (Path() 转换)。
  * R2-2-2  plot_arbitrary 状态着色: mask 与数据列同源, 含
            passive/lost 粒子的真实输出不再 ValueError。
  * R2-2-6  plot_phase_space 正电子正则散角: 用 canonical_signs
            的种类感知符号 (旧实现固定电子符号, 正电子散角为 0)。
"""

from pathlib import Path

import numpy as np
import pytest

from astra_tools.constants import C_LIGHT
from astra_tools.distribution import Distribution


# ---------------------------------------------------------------
# R2-2-1: phase_file_z_m
# ---------------------------------------------------------------

def test_phase_file_z_m_cm_convention():
    """3 位 cm 命名: Example.0150.001 -> 1.5 m."""
    from astra_tools.widgets.selectors import phase_file_z_m
    assert phase_file_z_m(Path("Example.0150.001")) == pytest.approx(1.5)


def test_phase_file_z_m_mm_convention():
    """4 位 mm 命名: Example.1500.001 -> 1.5 m."""
    from astra_tools.widgets.selectors import phase_file_z_m
    assert phase_file_z_m(Path("Example.1500.001")) == pytest.approx(1.5)


def test_phase_file_z_m_mid_dump():
    """中间 dump (螺线管中心场景): 1245 -> 1.245 m, 而非 ZSTOP 1.5 m."""
    from astra_tools.widgets.selectors import phase_file_z_m
    assert phase_file_z_m(Path("Example.1245.001")) == pytest.approx(1.245)


def test_phase_file_z_m_small_z():
    """小 z: 3 位按 cm 约定 (0050 -> 0.5 m, 0015 -> 0.15 m)."""
    from astra_tools.widgets.selectors import phase_file_z_m
    assert phase_file_z_m(Path("Example.0050.001")) == pytest.approx(0.5)
    assert phase_file_z_m(Path("Example.0015.001")) == pytest.approx(0.15)


def test_phase_file_z_m_rejects_non_phase_names():
    """非相位文件名 -> None (调用方显示告警, 不静默回退)。"""
    from astra_tools.widgets.selectors import phase_file_z_m
    assert phase_file_z_m(Path("bunch.ini")) is None
    assert phase_file_z_m(Path("astra.track.001")) is None
    assert phase_file_z_m(Path("README")) is None
    assert phase_file_z_m(Path("Example.Xemit.001")) is None


def test_phase_file_z_m_accepts_str():
    from astra_tools.widgets.selectors import phase_file_z_m
    assert phase_file_z_m("Example.0150.001") == pytest.approx(1.5)


def test_phase_label_accepts_str():
    """P3: _phase_label 接受 str (Path() 转换; 旧实现 f.name 崩溃)。"""
    from astra_tools.widgets.selectors import _phase_label
    assert _phase_label("Example.0150.001") == \
        "z = 1.5000 m  (Example.0150.001)"


def test_phase_label_regression_binary_cm_mm(tmp_path):
    """既有 Path 行为不回归: 二进制 (NUL 头) 文件走 cm/mm 文件名回退。"""
    from astra_tools.widgets.selectors import _phase_label
    f = tmp_path / "Example.0150.001"
    f.write_bytes(b"\x00\x01\x02\x03")
    assert _phase_label(f) == "z = 1.5000 m  (Example.0150.001)"
    f2 = tmp_path / "Example.1500.001"
    f2.write_bytes(b"\x00\x01\x02\x03")
    assert _phase_label(f2) == "z = 1.5000 m  (Example.1500.001)"


# ---------------------------------------------------------------
# R2-2-2: plot_arbitrary 状态着色 mask 与数据列同源
# ---------------------------------------------------------------

def _status_dist():
    """混合状态分布: active + passive + cathode + lost + lost deep."""
    n = 40
    y = np.linspace(-1e-3, 1e-3, n)
    status = np.full(n, 5, dtype=np.int32)
    status[10:20] = 0      # passive (探针)
    status[20:25] = -3     # cathode (未发射)
    status[25:30] = -10    # lost aperture
    status[30:35] = -30    # lost deep
    return Distribution.from_arrays(
        x=y, y=y, z=np.zeros(n), px=np.zeros(n), py=np.zeros(n),
        pz=np.full(n, 1e6), clock=np.zeros(n), charge=np.ones(n),
        status=status)


def test_plot_arbitrary_color_by_status_with_lost():
    """状态着色 + 含 passive/lost 粒子不抛 ValueError;
    曲线点数 = status>=-6 粒子数 (mask 与数据列同源)。"""
    import matplotlib.pyplot as plt
    from astra_tools.plot.arbitrary_phase_space import plot_arbitrary
    d = _status_dist()
    mask = d.status >= -6
    fig = plot_arbitrary(d, "x", "xp", color_by_status=True)
    try:
        pts = sum(len(c.get_offsets()) for c in fig.axes[0].collections)
        assert pts == int(mask.sum())
    finally:
        plt.close("all")


def test_plot_arbitrary_color_by_status_active_plus_passive():
    """active + passive (无 lost): 状态着色 mask 覆盖全部粒子;
    旧实现列来自 active 子集 -> 长度不匹配 ValueError。"""
    import matplotlib.pyplot as plt
    from astra_tools.plot.arbitrary_phase_space import plot_arbitrary
    n = 10
    status = np.array([5] * 5 + [0] * 5, dtype=np.int32)  # 5 active + 5 passive
    d = Distribution.from_arrays(
        x=np.zeros(n), y=np.zeros(n), z=np.zeros(n),
        px=np.zeros(n), py=np.zeros(n), pz=np.full(n, 1e6),
        clock=np.zeros(n), charge=np.ones(n), status=status)
    fig = plot_arbitrary(d, "x", "xp", color_by_status=True)
    try:
        pts = sum(len(c.get_offsets()) for c in fig.axes[0].collections)
        assert pts == n
    finally:
        plt.close("all")


def test_param_columns_mask_parameter():
    """param_columns 支持 mask 参数: 列长度 = mask 粒子数。"""
    from astra_tools.plot.arbitrary_phase_space import param_columns
    d = _status_dist()
    mask = d.status >= -6
    cols = param_columns(d, mask=mask)
    for name in ("x", "y", "xp", "yp", "E_kin", "dp/p"):
        assert len(cols[name][0]) == int(mask.sum()), name
    # 默认 (无 mask) 行为不变: active 粒子
    cols_active = param_columns(d)
    assert len(cols_active["x"][0]) == d.n_active


# ---------------------------------------------------------------
# R2-2-6: plot_phase_space 正电子正则散角符号
# ---------------------------------------------------------------

def test_plot_phase_space_positron_canonical_divergence():
    """正电子 (species=2) 用 canonical_signs 的符号: x' 非零。

    构造 px = -(c/2)·Bz·y (与电子符号项相反): 旧实现 (固定 +1)
    给出 p~x = 0 -> 全零散角; 正确符号 (s=-1) 给出 p~x = -c·Bz·y
    -> y=1 mm, Bz=0.1 T, p=1 MeV/c 时 x' = -30 mrad。
    """
    import matplotlib.pyplot as plt
    from astra_tools.plot.phase_space import plot_phase_space
    n = 50
    y = np.linspace(-1e-3, 1e-3, n)
    bz = 0.1
    px = -0.5 * C_LIGHT * bz * y
    d = Distribution.from_arrays(
        x=y, y=y, z=np.zeros(n), px=px, py=np.zeros(n),
        pz=np.full(n, 1e6), clock=np.zeros(n), charge=np.ones(n),
        index=np.full(n, 2, dtype=np.int32))
    fig = plot_phase_space(d, plane="x", bz_on_axis_T=bz)
    try:
        ys = np.concatenate(
            [c.get_offsets()[:, 1] for c in fig.axes[0].collections])
        assert np.std(ys) > 1e-6, \
            "正电子散角应为非零 (canonical_signs 种类感知符号)"
        # 量级检查: std(x') = c·Bz·std(y)/p_ref·1e3 ~ 17.3 mrad
        assert np.std(ys) > 10.0, "x' 量级应为 ~17 mrad, 实际 %r" % np.std(ys)
    finally:
        plt.close("all")


def test_plot_phase_space_electron_sign_unchanged():
    """电子 (species=1) 行为不回归: px = +(c/2)·Bz·y 构造 -> 散角非零。"""
    import matplotlib.pyplot as plt
    from astra_tools.plot.phase_space import plot_phase_space
    n = 50
    y = np.linspace(-1e-3, 1e-3, n)
    bz = 0.1
    px = +0.5 * C_LIGHT * bz * y
    d = Distribution.from_arrays(
        x=y, y=y, z=np.zeros(n), px=px, py=np.zeros(n),
        pz=np.full(n, 1e6), clock=np.zeros(n), charge=np.ones(n),
        index=np.full(n, 1, dtype=np.int32))
    fig = plot_phase_space(d, plane="x", bz_on_axis_T=bz)
    try:
        ys = np.concatenate(
            [c.get_offsets()[:, 1] for c in fig.axes[0].collections])
        assert np.std(ys) > 10.0, "电子散角不应回归为零"
    finally:
        plt.close("all")

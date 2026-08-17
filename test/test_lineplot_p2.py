"""lineplot P2 遗留项测试 (菜单 4 项 3/4/5/6, 菜单 1 项 15).

手册 4.13.6 (缩减发射度与 Eq. 4.4 K 项) 与 5.5.2 (探针轨迹柱坐标)。
"""

import numpy as np
import pytest

from astra_tools.distribution import Distribution
from astra_tools.plot.advanced_plots import (
    plot_correlated_emittance_contributions, plot_emittance_difference,
    plot_probe_trajectories, plot_reduced_longitudinal_emittance)


def _mk_dist(pz):
    n = len(pz)
    z = np.linspace(-1e-3, 1e-3, n)
    return Distribution.from_arrays(
        x=z, y=z * 0, z=z, px=z * 0, py=z * 0, pz=pz,
        clock=np.zeros(n), charge=np.ones(n))


def test_probe_trajectories_cylindrical():
    import matplotlib.pyplot as plt
    track = dict(seq=np.array([1, 1, 1, 2, 2, 2]),
                 z=np.array([0.0, 1.0, 2.0, 0.0, 1.0, 2.0]),
                 x=np.array([0.0, 1e-3, 2e-3, 0.0, 0.5e-3, 1e-3]),
                 y=np.array([0.0, 1e-3, 0.0, 0.0, 0.5e-3, 0.0]),
                 status=np.ones(6, int))
    fig = plot_probe_trajectories(track, mode="cylindrical")
    try:
        assert len(fig.axes) == 2
        assert "r [mm]" in fig.axes[0].get_ylabel()
        assert "x / y [mm]" in fig.axes[1].get_ylabel()
        # 第一探针 r(z=2) = hypot(2e-3, 0) = 2e-3 mm? -> 2.0 mm
        line = fig.axes[0].lines[0]
        np.testing.assert_allclose(line.get_ydata(), [0.0, np.sqrt(2), 2.0])
    finally:
        plt.close("all")


def test_emittance_difference_values():
    import matplotlib.pyplot as plt
    z = np.linspace(0, 1.5, 20)
    emit = dict(z=z, norm_emit_x=np.full(20, 2e-6), norm_emit_y=np.full(20, 3e-6))
    x2 = dict(z=z, eps_red_z=np.full(20, 1e-6), eps_red_zE=np.full(20, 1.2e-6))
    fig = plot_emittance_difference(emit, x2)
    try:
        labels = [ln.get_label() for ln in fig.axes[0].lines]
        assert any("x" in s for s in labels)
        for ln in fig.axes[0].lines:
            if "est" not in ln.get_label() and "y" in ln.get_label():
                np.testing.assert_allclose(ln.get_ydata(), 2.0, rtol=1e-6)
        # x: 2e-6 - 1e-6 = 1e-6 -> 1.0 [pi mm mrad]
        xline = [ln for ln in fig.axes[0].lines if "x" in ln.get_label()
                 and "est" not in ln.get_label()][0]
        np.testing.assert_allclose(xline.get_ydata(), 1.0, rtol=1e-6)
    finally:
        plt.close("all")


def test_correlated_contributions_values():
    import matplotlib.pyplot as plt
    z = np.linspace(0, 1.5, 20)
    x2 = dict(z=z, K2z=np.full(20, 2.0), K3z=np.full(20, 0.5),
              K2E=np.full(20, 1.0), K3E=np.full(20, 0.25),
              eps_red_z=np.zeros(20), eps_red_zE=np.zeros(20))
    fig = plot_correlated_emittance_contributions(x2)
    try:
        by_label = {ln.get_label(): ln.get_ydata() for ln in fig.axes[0].lines}
        # K 项文件数值 rad.m, 显示 *1e6 -> [pi mm mrad]
        assert np.allclose(by_label.get("$K_{2,Z}$ [x]", [np.nan]), 2e6)
        assert np.allclose(by_label.get("$K_{3,E}$ [x]", [np.nan]), 0.25e6)
    finally:
        plt.close("all")


def test_reduced_longitudinal_uncorrelated():
    """无 z-pz 相关: 去相关后发射度几乎不变."""
    import matplotlib.pyplot as plt
    rng = np.random.default_rng(0)
    n = 3000
    z = np.linspace(-1e-3, 1e-3, n)
    pz = 1e6 + rng.normal(0, 1e3, n)
    fig = plot_reduced_longitudinal_emittance(_mk_dist(pz))
    try:
        txt = fig.axes[0].texts[0].get_text()
        import re
        vals = [float(v) for v in re.findall(r"= ([\d.eE+-]+)", txt)]
        assert len(vals) == 2
        assert vals[1] == pytest.approx(vals[0], rel=0.2)
    finally:
        plt.close("all")


def test_reduced_longitudinal_strong_corr():
    """强二次 z-pz 相关: 去相关 (3rd 阶含 z²) 后发射度显著减小.

    几何发射度公式自动消除线性相关, 因此用二次相关构造测试 (手册
    4.13.6: RF 场引入 C2 二次相关)。
    """
    import matplotlib.pyplot as plt
    rng = np.random.default_rng(1)
    n = 3000
    z = np.linspace(-1e-3, 1e-3, n)
    pz = 1e6 + 2e9 * z ** 2 + rng.normal(0, 1e2, n)
    fig = plot_reduced_longitudinal_emittance(_mk_dist(pz))
    try:
        txt = fig.axes[0].texts[0].get_text()
        import re
        vals = [float(v) for v in re.findall(r"= ([\d.eE+-]+)", txt)]
        assert len(vals) == 2
        assert vals[1] < 0.3 * vals[0]
    finally:
        plt.close("all")


def test_emit_diff_y_est_without_yemit():
    """无 Yemit2 时 y 平面用 x 近似并标注 est.."""
    import matplotlib.pyplot as plt
    z = np.linspace(0, 1.5, 20)
    emit = dict(z=z, norm_emit_x=np.full(20, 2e-6), norm_emit_y=np.full(20, 4e-6))
    x2 = dict(z=z, eps_red_z=np.full(20, 1e-6), eps_red_zE=np.full(20, 1.2e-6))
    fig = plot_emittance_difference(emit, x2)
    try:
        labels = [ln.get_label() for ln in fig.axes[0].lines]
        assert any("est" in s for s in labels)
    finally:
        plt.close("all")

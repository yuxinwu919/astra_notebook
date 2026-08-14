"""绘图正确性测试 (第 4 层): 每个绘图函数的单位标签/图例/数据完整性.

历史上单位与物理量是绘图代码最容易出错的地方, 本测试把关键
轴标签与单位字符串固化为断言, 并确保真实数据上无 NaN 密度。
"""

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.io import read_distribution
from astra_tools.io.astra_emit import read_emit_files, read_ref_file
from astra_tools.io.field_map import read_cavity_field, read_solenoid_field
from astra_tools.analysis.slices import compute_slice_analysis
from astra_tools.analysis.bff import compute_bff
from astra_tools.plot.style import set_style
from astra_tools.plot.phase_space import plot_phase_space
from astra_tools.plot.distributions import plot_distributions, plot_energy_distribution
from astra_tools.plot.overview import plot_overview, plot_transverse_profile
from astra_tools.plot.emit_plots import (
    plot_emittance_evolution, plot_energy_evolution, plot_emit_dashboard,
    plot_lineplot_overview, plot_ref_trajectory,
)
from astra_tools.plot.slice_plots import plot_slice_dashboard
from astra_tools.plot.bff_plots import plot_bff
from astra_tools.plot.field_plots import plot_cavity_field, plot_solenoid_field

set_style()

DIST = read_distribution(PROJECT_ROOT / "examples/Manual_Example/Example.0150.001")
EMIT = read_emit_files(str(PROJECT_ROOT / "examples/Manual_Example/Example"))
REF = read_ref_file(str(PROJECT_ROOT / "examples/Manual_Example/Example"))
SA = compute_slice_analysis(DIST, n_slices=10)
BFF = compute_bff(DIST.filter_active().z, DIST.filter_active().charge)
CAV = read_cavity_field(PROJECT_ROOT / "examples/Manual_Example/3_cell_L-Band.dat")
SOL = read_solenoid_field(PROJECT_ROOT / "examples/Manual_Example/Solenoid.dat").scaled(0.35)


def _labels(fig):
    out = []
    for ax in fig.axes:
        out += [ax.get_xlabel(), ax.get_ylabel()]
    return out


class TestPhaseSpace:
    def test_transverse_units(self):
        fig = plot_phase_space(DIST, plane="x", show_ellipse=True)
        labs = _labels(fig)
        assert "x [mm]" in labs
        assert "x' [mrad]" in labs
        ax = fig.axes[0]
        assert ax.get_legend() is not None, "RMS 椭圆必须有图例"
        z = ax.collections[0].get_array()
        assert np.all(np.isfinite(z)), "密度数据不得含 NaN"
        plt.close(fig)

    def test_longitudinal_units(self):
        fig = plot_phase_space(DIST, plane="z")
        labs = _labels(fig)
        assert any("dp/p [%]" in l for l in labs)
        assert any("z [mm]" in l for l in labs)
        plt.close(fig)

    def test_outlier_clip_annotation(self):
        n = 2010
        rng = np.random.default_rng(0)
        x = np.concatenate([rng.normal(0, 1e-3, 2000),
                            np.full(10, 0.5)])  # 10 个 0.5 m 极端离群点 (0.5%)
        from astra_tools.distribution import Distribution
        d = Distribution.from_arrays(
            x, np.zeros(n), np.zeros(n), np.zeros(n), np.zeros(n),
            np.full(n, 1e6), np.zeros(n), np.full(n, 1e-3),
            status=np.full(n, 5), ref_momentum_eVc=1e6)
        fig = plot_phase_space(d, plane="x")
        texts = [t.get_text() for t in fig.axes[0].texts]
        assert any("range clip" in t for t in texts), "离群点裁剪需有标注"
        xl = fig.axes[0].get_xlim()
        # 轴单位为 mm: 主体 ±3 mm 可见, 500 mm 离群点被裁掉
        assert abs(xl[0]) <= 5.0 and xl[1] <= 5.0, "主体分布必须可见(毫米级), 实际 xlim=%s" % (xl,)
        plt.close(fig)


class TestEmitPlots:
    def test_emittance_unit_is_pi_mm_mrad(self):
        fig = plot_emittance_evolution(EMIT)
        yl = fig.axes[0].get_ylabel()
        assert "$\\pi$ mm mrad" in yl, "发射度轴必须为 pi mm mrad (与 ASTRA 同步)"
        plt.close(fig)

    def test_energy_evolution(self):
        fig = plot_energy_evolution(EMIT)
        labs = _labels(fig)
        assert any("[MeV]" in l for l in labs)
        assert any("[keV]" in l for l in labs)
        plt.close(fig)

    def test_overview_panels(self):
        fig = plot_lineplot_overview(EMIT)
        assert len(fig.axes) >= 9  # 9 主面板 (+ twin 轴)
        plt.close(fig)

    def test_ref_trajectory(self):
        fig = plot_ref_trajectory(REF)
        labs = _labels(fig)
        assert any("[MeV]" in l for l in labs)
        plt.close(fig)


class TestOtherPlots:
    def test_distributions(self):
        fig = plot_distributions(DIST)
        assert len(fig.axes) == 3
        for ax in fig.axes:
            assert ax.get_legend() is not None
        plt.close(fig)

    def test_overview_grid(self):
        fig, axes = plot_overview(DIST)
        assert axes.shape == (3, 2)
        plt.close(fig)

    def test_slice_dashboard(self):
        fig = plot_slice_dashboard(SA)
        assert len(fig.axes) >= 4  # 4 主面板 (+ twin 轴)
        plt.close(fig)

    def test_bff(self):
        fig = plot_bff(BFF)
        assert fig.axes[0].get_xlabel() == "k [1/m]"
        plt.close(fig)

    def test_cavity_field_map(self):
        fig = plot_cavity_field(CAV, maxE_MVpm=40)
        labs = _labels(fig)
        assert any("[MV/m]" in l for l in labs)
        plt.close(fig)

    def test_solenoid_field_map(self):
        fig = plot_solenoid_field(SOL)
        labs = _labels(fig)
        assert any("[T]" in l for l in labs)
        plt.close(fig)

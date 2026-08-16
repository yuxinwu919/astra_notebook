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

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.io import read_distribution
from astra_tools.io.astra_emit import read_emit_files, read_ref_file
from astra_tools.io.field_map import read_cavity_field, read_solenoid_field
from astra_tools.analysis.slices import compute_slice_analysis
from astra_tools.analysis.bff import compute_bff
from astra_tools.plot.style import set_style
from astra_tools.plot.phase_space import plot_phase_space
from astra_tools.plot.distributions import plot_distributions
from astra_tools.plot.overview import plot_overview
from astra_tools.plot.emit_plots import (
    plot_emittance_evolution, plot_energy_evolution,
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
        fig = plot_phase_space(DIST, plane="x")
        labs = _labels(fig)
        assert "x [mm]" in labs
        assert "x' [mrad]" in labs
        ax = fig.axes[0]
        assert len(ax.collections) >= 1, "散点集合必须存在"
        z = ax.collections[0].get_offsets()
        assert np.all(np.isfinite(z)), "散点数据不得含 NaN"
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


class TestAdvancedPlots:
    """lineplot 菜单 2/3/4 与 postpro 扩展绘图."""

    def test_beta_alpha_units(self):
        from astra_tools.plot.advanced_plots import plot_beta_alpha
        fig = plot_beta_alpha(EMIT)
        labs = _labels(fig)
        assert any("beta function [m]" in l for l in labs)
        assert fig.axes[0].get_legend() is not None
        plt.close(fig)

    def test_phase_advance(self):
        from astra_tools.plot.advanced_plots import plot_phase_advance
        fig = plot_phase_advance(EMIT)
        assert any("phase advance [rad]" in l for l in _labels(fig))
        plt.close(fig)

    def test_coherence_length_units(self):
        from astra_tools.plot.advanced_plots import plot_coherence_length
        fig = plot_coherence_length(EMIT)
        assert any("coherence length [m]" in l for l in _labels(fig))
        plt.close(fig)

    def test_slice_mismatch_ge_one(self):
        """失配参数物理下界: zeta >= 1 (手册 5.6.3)."""
        from astra_tools.plot.advanced_plots import slice_mismatch
        _, zx, zy = slice_mismatch(DIST, n_slices=8)
        assert np.all(np.nanmin(zx) >= 0.999)
        assert np.all(np.nanmin(zy) >= 0.999)

    def test_3d_map_slices(self):
        from astra_tools.plot.advanced_plots import plot_3d_map_slices
        p = PROJECT_ROOT / "examples/Cavity_Example/3D_test.ex"
        fig = plot_3d_map_slices(p, axis="z", n_slices=2, unit="V/m")
        assert len(fig.axes) >= 2
        plt.close(fig)

    def test_3d_map_reader(self):
        from astra_tools.io.field_map import read_3d_field_map
        x, y, z, f = read_3d_field_map(
            PROJECT_ROOT / "examples/Cavity_Example/3D_test.ex")
        assert f.shape == (11, 11, 340)
        assert np.all(np.isfinite(f))


class TestCutsAndMisc:
    def test_cut_window(self):
        from astra_tools.analysis.cuts import cut_distribution
        d2, mask = cut_distribution(DIST, x_range=(-5e-4, 5e-4))
        assert mask.sum() > 0
        assert d2.n_active == DIST.n_active - mask.sum()
        assert np.all(d2.status[mask] == -31)

    def test_rotate_phase_space(self):
        from astra_tools.analysis.cuts import rotate_phase_space
        d3 = rotate_phase_space(DIST, 90)
        assert d3.n_particle == DIST.n_particle
        # 旋转 90 度: x -> y
        assert np.allclose(d3.x[DIST.active], DIST.y[DIST.active], atol=1e-12)

    def test_z_plot(self):
        from astra_tools.plot.advanced_plots import plot_z_plot
        fig = plot_z_plot(DIST)
        ax = fig.axes[0]
        assert ax.get_legend() is not None
        assert ax.get_ylabel() == "z [mm]"
        plt.close(fig)

    def test_slice_ellipses_3d(self):
        from astra_tools.plot.advanced_plots import plot_slice_ellipses_3d
        fig = plot_slice_ellipses_3d(DIST, n_slices=6)
        assert fig.axes[0].name == "3d"
        plt.close(fig)

    def test_cathode_contour(self):
        from astra_tools.plot.advanced_plots import plot_curved_cathode_contour
        p = PROJECT_ROOT / "examples/Curved_Cathode_Example/Contour.dat"
        fig = plot_curved_cathode_contour(p)
        assert any("[mm]" in l for l in _labels(fig))
        plt.close(fig)

    def test_font_policy_times_new_roman(self):
        """绘图字体统一 Times New Roman (用户决定 2026-08), 回退链兜底."""
        fam = plt.rcParams["font.family"]
        assert fam[0] == "Times New Roman"
        assert "STIXGeneral" in fam and "DejaVu Sans" in fam
        assert plt.rcParams["mathtext.fontset"] == "stix"
        # 回退链里绝不能残留本机不存在的字体: 否则每个文字对象都会
        # 为每个缺失家族发一条 findfont 警告 (数千条 IOPub 消息,
        # Jupyter 单元会被拖到几分钟不结束)。
        from matplotlib import font_manager
        installed = {f.name for f in font_manager.fontManager.ttflist}
        assert all(f in installed for f in fam)

    def test_font_chain_covers_every_used_glyph(self):
        """包内所有字符串字面量用到的每个非 ASCII 字符, 回退链都必须覆盖.

        "字体显示不全"的根治验证: 逐个字形查 charmap, 任何字符无字体
        覆盖即为方框来源。mathtext 字符 (\beta 等) 由 STIX 字库渲染,
        不在此检查范围。"""
        import ast
        import pathlib
        from matplotlib import font_manager
        from matplotlib.ft2font import FT2Font
        fam = plt.rcParams["font.family"]
        installed = {f.name: f.fname for f in font_manager.fontManager.ttflist}
        cmaps = {}
        for name in fam:
            if name in installed:
                cmaps[name] = set(FT2Font(installed[name]).get_charmap().keys())
        used = set()
        for f in pathlib.Path(PROJECT_ROOT, "astra_tools").rglob("*.py"):
            try:
                tree = ast.parse(f.read_text())
            except SyntaxError:
                continue
            for node in ast.walk(tree):
                if isinstance(node, ast.Constant) and isinstance(node.value, str):
                    used |= {c for c in node.value if ord(c) > 127}
        uncovered = [c for c in sorted(used)
                     if not any(ord(c) in cm for cm in cmaps.values())]
        assert not uncovered, "无字体覆盖的字形: %s" % "".join(uncovered)

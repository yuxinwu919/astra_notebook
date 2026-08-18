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
        fig = plot_beta_alpha(EMIT, ref=REF)
        labs = _labels(fig)
        assert any("beta function [m]" in l for l in labs)
        assert fig.axes[0].get_legend() is not None
        # 批 2: beta 用几何发射度, 数值按 sigma^2/eps_geom 给出
        plt.close(fig)

    def test_phase_advance(self):
        from astra_tools.plot.advanced_plots import plot_phase_advance
        fig = plot_phase_advance(EMIT, ref=REF)
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
        from astra_tools.plot.field_plots import plot_3d_field_map
        p = PROJECT_ROOT / "examples/Cavity_Example/3D_test.ex"
        fig = plot_3d_field_map(p, view="slices", component="z", n_slices=2)
        assert len(fig.axes) >= 2
        plt.close(fig)

    def test_3d_map_slices_all_zero_field_warns(self, tmp_path):
        """全零场切面应告警 (常见原因: 选错分量文件), 而非静默空白图。"""
        from astra_tools.plot.field_plots import plot_3d_field_map
        p = tmp_path / "zero_field.ex"
        p.write_text(
            "3 0.0 1.0 2.0\n"
            "3 0.0 1.0 2.0\n"
            "3 0.0 1.0 2.0\n"
            + " ".join(["0.0"] * 27) + "\n")
        with pytest.warns(UserWarning, match="全为零"):
            fig = plot_3d_field_map(p, view="slices", n_slices=3)
        assert len(fig.axes) >= 3
        plt.close(fig)

    def test_3d_map_reader(self):
        from astra_tools.io.field_map import read_3d_field_map
        x, y, z, f = read_3d_field_map(
            PROJECT_ROOT / "examples/Cavity_Example/3D_test.ex")
        assert f.shape == (11, 11, 340)
        assert np.all(np.isfinite(f))


class Test3DFieldMapViews:
    """FieldMap3D 数据层与 plot_3d_field_map 四种视图 (fieldplot 菜单 2)."""

    DIPOLE = PROJECT_ROOT / "examples/90deg_bend_Example/3D_Dipole"

    def test_read_components_units_and_magnitude(self):
        from astra_tools.io import read_3d_field_map_components
        f = read_3d_field_map_components(self.DIPOLE)
        assert f.unit == "T" and f.quantity == "B"
        assert f.fx.shape == (40, 3, 46)
        # Bx=Bz=0 -> 模长 = |By|
        assert np.allclose(f.magnitude, np.abs(f.fy))
        # 显式后缀仍可推导单位
        assert read_3d_field_map_components(
            str(self.DIPOLE) + ".by").unit == "T"
        # 主名无后缀且电/磁分量并存时优先电场
        cav = read_3d_field_map_components(
            PROJECT_ROOT / "examples/Cavity_Example/3D_test")
        assert cav.unit == "V/m" and cav.quantity == "E"
        # 电场分量不被磁场文件覆盖 (3D_test 的 bx/by/bz 是噪声, 全零/1e-8)
        assert np.max(np.abs(cav.fx)) > 1.0
        assert np.max(np.abs(cav.fz)) > 1.0
        # 显式 .ez 后缀与主名电场一致
        cav_ez = read_3d_field_map_components(
            str(PROJECT_ROOT / "examples/Cavity_Example/3D_test") + ".ez")
        assert cav_ez.unit == "V/m"
        assert np.allclose(cav_ez.fz, cav.fz)

    def test_read_components_missing_and_mismatch(self, tmp_path):
        from astra_tools.io import read_3d_field_map_components

        def write(name, grid, val):
            p = tmp_path / name
            p.write_text(
                "3 %s\n" % " ".join(str(g) for g in grid) * 3
                + " ".join([str(val)] * 27) + "\n")

        write("map.bx", (0.0, 0.5, 1.0), 1.0)
        write("map.by", (0.0, 0.5, 1.0), 2.0)
        write("map.bz", (0.0, 0.5, 1.0), 0.0)
        f = read_3d_field_map_components(tmp_path / "map")
        assert f.unit == "T"
        assert np.allclose(f.magnitude, np.sqrt(5.0))
        # 缺失分量按全零
        (tmp_path / "map.bz").unlink()
        f2 = read_3d_field_map_components(tmp_path / "map")
        assert np.all(f2.fz == 0.0)
        # 网格不一致必须报错
        write("mx.bx", (0.0, 0.5, 1.0), 1.0)
        write("mx.by", (0.0, 1.0, 2.0), 2.0)
        with pytest.raises(ValueError, match="网格不一致"):
            read_3d_field_map_components(tmp_path / "mx")

    def test_vector_slices_labels(self):
        """矢量剖面: auto 轴 (y) + mm 标签 + |B| [T] 色条。"""
        from astra_tools.plot.field_plots import plot_3d_field_map
        fig = plot_3d_field_map(self.DIPOLE, view="slices",
                                n_slices=2)
        labs = _labels(fig)
        assert any("x [mm]" in l for l in labs)
        assert any("z [mm]" in l for l in labs)
        assert "|B| [T]" in (fig.axes[-1].get_ylabel() or "")
        plt.close(fig)

    def test_plane_labels_and_z_horizontal(self):
        """plane 语义 + xz/yz 平面中 z 在横轴约定。"""
        from astra_tools.plot.field_plots import plot_3d_field_map
        fig = plot_3d_field_map(self.DIPOLE, view="slices",
                                plane="xy", n_slices=2)
        assert fig.axes[0].get_xlabel() == "x [mm]"
        assert fig.axes[0].get_ylabel() == "y [mm]"
        plt.close(fig)
        fig = plot_3d_field_map(self.DIPOLE, view="slices",
                                plane="xz", n_slices=2)
        assert fig.axes[0].get_xlabel() == "z [mm]"     # z 在横轴
        assert fig.axes[0].get_ylabel() == "x [mm]"
        plt.close(fig)
        fig = plot_3d_field_map(self.DIPOLE, view="slices",
                                plane="yz", n_slices=2)
        assert fig.axes[0].get_xlabel() == "z [mm]"
        assert fig.axes[0].get_ylabel() == "y [mm]"
        plt.close(fig)

    def test_single_plane_selection(self):
        """position/index 指定单个剖面 (单面板)。"""
        from astra_tools.plot.field_plots import plot_3d_field_map
        fig = plot_3d_field_map(self.DIPOLE, view="slices",
                                plane="xz", position=0.0)
        assert len(fig.axes) == 2
        assert fig.axes[0].get_title() == "y = 0 mm"
        plt.close(fig)
        fig = plot_3d_field_map(self.DIPOLE, view="slices", kind="contour",
                                plane="xy", index=22)
        assert fig.axes[0].get_title() == "z = -5 mm"
        plt.close(fig)
        with pytest.raises(ValueError, match="越界"):
            plot_3d_field_map(self.DIPOLE, view="slices",
                              plane="xy", index=999)

    def test_interactive_slices_widget(self):
        """滑块组件: 构造 + 滑块范围与固定轴网格一致 (auto -> xz/y)。"""
        from astra_tools.widgets import interact_3d_field_slices
        wb = interact_3d_field_slices(self.DIPOLE, auto_render=False)
        slider = wb.children[2]
        assert abs(slider.min - (-3.0)) < 1e-9
        assert abs(slider.max - 3.0) < 1e-9
        assert wb.children[0].value == "heatmap"
        assert wb.children[1].value == "xz"

    @staticmethod
    def _quiver_collection(fig):
        from matplotlib.quiver import Quiver
        return next(c for c in fig.axes[0].collections
                    if isinstance(c, Quiver))

    def test_quiver_2d_semantics(self):
        """纯 2D quiver: plane='xy' -> xy 平面, 箭头 (U,V) = (fx, fy)."""
        from astra_tools.io import read_3d_field_map_components
        from astra_tools.plot.field_plots import plot_3d_field_quiver
        f = read_3d_field_map_components(self.DIPOLE)
        fig = plot_3d_field_quiver(f, plane="xy", index=22)
        ax = fig.axes[0]
        assert ax.get_xlabel() == "x [mm]"
        assert ax.get_ylabel() == "y [mm]"
        assert ax.get_title() == "z = -5 mm"
        q = self._quiver_collection(fig)
        assert np.asarray(q.U).shape == (3 * 40,)    # matplotlib 展平存储
        assert np.allclose(np.asarray(q.U), f.fx[:, :, 22].T.ravel())
        assert np.allclose(np.asarray(q.V), f.fy[:, :, 22].T.ravel())
        plt.close(fig)

    def test_quiver_axis_semantics(self):
        """plane='xz' (z 横轴, U,V)=(fz,fx); plane='yz' (fz,fy)."""
        from astra_tools.io.field_map import FieldMap3D
        from astra_tools.plot.field_plots import plot_3d_field_quiver
        # 合成三分量非零场: fx=x, fy=y, fz=z (每个剖面面内分量明确)
        x = np.linspace(-1, 1, 5)
        y = np.linspace(-2, 2, 4)
        z = np.linspace(-3, 3, 6)
        XX, YY, ZZ = np.meshgrid(x, y, z, indexing="ij")
        f = FieldMap3D(x=x, y=y, z=z, fx=XX, fy=YY, fz=ZZ,
                       unit="V/m", quantity="E")
        fig = plot_3d_field_quiver(f, plane="xz", index=1)
        ax = fig.axes[0]
        assert ax.get_xlabel() == "z [mm]"
        assert ax.get_ylabel() == "x [mm]"
        q = self._quiver_collection(fig)
        assert np.allclose(np.asarray(q.U), f.fz[:, 1, :].ravel())
        assert np.allclose(np.asarray(q.V), f.fx[:, 1, :].ravel())
        plt.close(fig)
        fig = plot_3d_field_quiver(f, plane="yz", index=0)
        ax = fig.axes[0]
        assert ax.get_xlabel() == "z [mm]"
        assert ax.get_ylabel() == "y [mm]"
        q = self._quiver_collection(fig)
        assert np.allclose(np.asarray(q.U), f.fz[0, :, :].ravel())
        assert np.allclose(np.asarray(q.V), f.fy[0, :, :].ravel())
        plt.close(fig)
        # 抽稀: max_arrows 上限 -> 等步长且与源数据一致
        x2 = np.linspace(0, 1, 20)
        y2 = np.linspace(0, 2, 30)
        z2 = np.linspace(0, 3, 40)
        X2, Y2, Z2 = np.meshgrid(x2, y2, z2, indexing="ij")
        f2 = FieldMap3D(x=x2, y=y2, z=z2, fx=X2, fy=Y2, fz=Z2,
                        unit="V/m", quantity="E")
        fig = plot_3d_field_quiver(f2, plane="xz", index=15, max_arrows=100)
        step = int(np.ceil(np.sqrt(f2.fz[:, 15, :].size / 100)))
        q = self._quiver_collection(fig)
        assert np.asarray(q.U).size <= 100
        assert np.allclose(np.asarray(q.U),
                           f2.fz[:, 15, :][::step, ::step].ravel())
        plt.close(fig)

    def test_quiver_color_and_validation(self):
        """color_by='magnitude' 加色条; 未知 plane/color_by 报错; 默认中间层. """
        from astra_tools.io import read_3d_field_map_components
        from astra_tools.plot.field_plots import plot_3d_field_quiver
        f = read_3d_field_map_components(self.DIPOLE)
        fig = plot_3d_field_quiver(f, color_by="magnitude")
        assert len(fig.axes) == 2                     # 图 + 色条
        assert "|B| [T]" in (fig.axes[1].get_ylabel() or "")
        assert fig.axes[0].get_title() == "z = %.4g mm" % (
            f.z[len(f.z) // 2] * 1e3)
        plt.close(fig)
        with pytest.raises(ValueError, match="plane"):
            plot_3d_field_quiver(f, plane="w")
        with pytest.raises(ValueError, match="color_by"):
            plot_3d_field_quiver(f, color_by="bogus")
        # 直接传路径与传 FieldMap3D 等价 (统一入口)
        fig2 = plot_3d_field_quiver(self.DIPOLE, plane="xy", index=22)
        assert np.allclose(
            np.asarray(self._quiver_collection(fig2).U),
            f.fx[:, :, 22].T.ravel())
        plt.close(fig2)

    def test_quiver_zero_plane_warns(self, tmp_path):
        """面内分量全零剖面: 告警而非崩溃 (如 3D_Dipole 的 xz 平面)."""
        import warnings
        from astra_tools.plot.field_plots import plot_3d_field_quiver
        # 全零场: 任何平面面内分量都为 0
        p = tmp_path / "zero.ex"
        p.write_text("3 0.0 0.5 1.0\n" * 3 + " ".join(["0.0"] * 27) + "\n")
        with pytest.warns(UserWarning, match="全为零"):
            fig = plot_3d_field_quiver(p, plane="xy")
        assert fig.axes[0].get_xlabel() == "x [mm]"
        assert len(fig.axes[0].collections) == 0      # 无 quiver
        plt.close(fig)
        # 但 Dipole 的 xy 平面 (U=fx=0, V=fy=By) 非零, 不告警
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            fig = plot_3d_field_quiver(self.DIPOLE, plane="xy", index=22)
        plt.close(fig)

    def test_interactive_quiver_widget(self):
        """quiver 滑块组件: 构造 + 固定平面按钮 + 滑块范围 = 该轴网格. """
        from astra_tools.io import read_3d_field_map_components
        from astra_tools.widgets import interact_3d_field_quiver
        f = read_3d_field_map_components(self.DIPOLE)
        wb = interact_3d_field_quiver(self.DIPOLE, plane="xy", auto_render=False)
        assert wb.children[0].value == "xy"
        assert wb.children[0].options == ("xy", "xz", "yz")
        slider = wb.children[1]
        assert abs(slider.min - float(f.z.min() * 1e3)) < 1e-9
        assert abs(slider.max - float(f.z.max() * 1e3)) < 1e-9
        assert slider.step > 0

    def test_contour_and_scalar_labels(self):
        from astra_tools.plot.field_plots import plot_3d_field_map
        fig = plot_3d_field_map(self.DIPOLE, view="slices", kind="contour",
                                n_slices=2)
        assert any("mm" in l for l in _labels(fig))
        assert "|B| [T]" in (fig.axes[-1].get_ylabel() or "")
        plt.close(fig)
        fig = plot_3d_field_map(self.DIPOLE, view="slices",
                                component="y", n_slices=2)
        assert "$B_y$ [T]" in (fig.axes[-1].get_ylabel() or "")
        plt.close(fig)

    def test_stack3d_renders(self):
        from astra_tools.plot.field_plots import plot_3d_field_map
        fig = plot_3d_field_map(self.DIPOLE, view="stack3d", n_slices=3)
        assert any(a.name == "3d" for a in fig.axes)
        assert any("x [mm]" in l for l in _labels(fig))
        plt.close(fig)

    def test_stack3d_aspect_modes(self):
        """3D 盒子比例: auto 模式保证最短轴有可读性下限 (不压扁)。"""
        from astra_tools.io import read_3d_field_map_components
        from astra_tools.plot.field_plots import _box_aspect
        f = read_3d_field_map_components(self.DIPOLE)
        phys = _box_aspect(f, "physical")
        assert abs(phys[1] / phys[2] - 6.0 / 450.0) < 1e-6  # y 被压扁
        assert _box_aspect(f, "equal") == (1.0, 1.0, 1.0)
        auto = _box_aspect(f, "auto")
        assert auto[1] >= 0.25 * auto[2]                      # 可读性下限
        assert _box_aspect(f, "grid") == (40.0, 3.0, 46.0)
        assert _box_aspect(f, (2.0, 1.0, 3.0)) == (2.0, 1.0, 3.0)
        with pytest.raises(ValueError, match="aspect"):
            _box_aspect(f, "bogus")

    def test_dispatcher_validation(self):
        from astra_tools.plot.field_plots import plot_3d_field_map
        with pytest.raises(ValueError, match="未知 view"):
            plot_3d_field_map(self.DIPOLE, view="nope")
        # 旧 view 字符串 (scalar_slices 等) 现在一律拒绝为未知 view
        with pytest.raises(ValueError, match="未知 view"):
            plot_3d_field_map(self.DIPOLE, view="scalar_slices")

    def test_zero_field_warns(self, tmp_path):
        """全零场: 矢量剖面应告警而非静默空白图。"""
        from astra_tools.plot.field_plots import plot_3d_field_map
        p = tmp_path / "zero.ex"
        p.write_text("3 0.0 0.5 1.0\n" * 3 + " ".join(["0.0"] * 27) + "\n")
        with pytest.warns(UserWarning, match="全为零"):
            fig = plot_3d_field_map(p, view="slices")
        plt.close(fig)


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
        line = fig.axes[0].lines[0]
        x = np.asarray(line.get_xdata(), dtype=float)
        y = np.asarray(line.get_ydata(), dtype=float)
        # 手册 4.4.5: 前两列 (z, R) 是坐标, 后两列是切向单位矢量。
        # 若误用切向分量 (范围 -1..1) 会得到 ~1000 mm 的假轮廓, 这里限定物理量级。
        assert -50 < x.min() and x.max() < 50, \
            "z 超出合理阴极尺寸 (误用了切向分量列?): %r" % ((x.min(), x.max()),)
        assert 0 <= y.min() and y.max() < 50, \
            "r 超出合理阴极尺寸 (误用了切向分量列?): %r" % ((y.min(), y.max()),)
        plt.close(fig)

    def test_trace_emittance_includes_longitudinal(self):
        """TRemit 含 eps_tr_z 时显示纵向面板 (菜单 4 项 8)."""
        from astra_tools.plot.advanced_plots import plot_trace_emittance
        z = np.linspace(0, 1.5, 20)
        tr = dict(z=z, eps_tr_x=np.full(20, 1e-6),
                  eps_tr_y=np.full(20, 1e-6), eps_tr_z=np.full(20, 1.2e-6))
        fig = plot_trace_emittance(tr)
        assert any("um" in l for l in _labels(fig)), "缺纵向 Trace 单位标签"
        plt.close(fig)

    def test_cr_emit_beam_size_panel(self):
        """Cr_emit 含 x_rms/y_rms 时画束斑面板 (菜单 4 项 13)."""
        from astra_tools.plot.advanced_plots import plot_cr_emit
        z = np.linspace(0, 1.5, 20)
        cr = dict(z=z, eps_x=np.full(20, 1e-6), eps_y=np.full(20, 1e-6),
                  q_rest=np.linspace(1, 0.8, 20), q_cross=np.linspace(0, 0.2, 20),
                  x_rms=np.linspace(1e-3, 2e-3, 20),
                  y_rms=np.linspace(1e-3, 1.5e-3, 20))
        fig = plot_cr_emit(cr)
        assert len(fig.axes) >= 2, "应含发射度 + 束斑两面板"
        assert any("beam size" in l.lower() for l in _labels(fig)), \
            "缺排除 cross-over 束斑面板"
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

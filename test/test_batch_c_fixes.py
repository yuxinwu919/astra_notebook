"""Batch C regression tests: 手册第5章末梢展示 + 真实 PScan/Scan 交叉验证.

PScan / Scan 部分是对真实 ASTRA 运行的独立交叉验证 (golden 文件由
本地 ASTRA 真实生成, 测试本身不需要可执行文件):

  * PScan: 单腔余弦律拟合残差、峰谷对称性、压缩因子、以及与同一
    运行 ref.001 (参考粒子) 末态能量的吻合 (<0.5%; 差异来自扫描
    粒子与分布参考粒子之间约 2° 的注入相位偏移)
  * Scan: FOM(1)=归一化水平发射度、FOM(2)=束长、FOM(3)=平均能量,
    与同一算例的 golden Xemit/Zemit 逐列对照 (MaxB=0.35 T 时)
"""

from pathlib import Path
import sys

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

DATA = PROJECT_ROOT / "examples" / "Manual_Example"
CAV = PROJECT_ROOT / "examples" / "Cavity_Example" / "golden"

from astra_tools.io import read_distribution
from astra_tools.io.astra_misc import read_pscan, read_scan
from astra_tools.io.astra_emit import read_emit_files, read_ref_file
from astra_tools.analysis.statistics import compute_statistics
from astra_tools.analysis.core import compute_central_charge_fraction_curves
import matplotlib
matplotlib.use("Agg")


# ---------------------------------------------------------------
# PScan 真实数据交叉验证
# ---------------------------------------------------------------
@pytest.fixture(scope="module")
def pscan():
    return read_pscan(CAV / "astra.PScan.001")


def test_pscan_single_cavity_cosine_law(pscan):
    """E(phi) = E0 + A cos(phi - phi0): 拟合残差 < 0.1% A."""
    ph = np.radians(pscan["phase_deg"])
    A = np.column_stack([np.ones_like(ph), np.cos(ph), np.sin(ph)])
    c, *_ = np.linalg.lstsq(A, pscan["E_kin_eV"], rcond=None)
    res = pscan["E_kin_eV"] - A @ c
    amp = np.hypot(c[1], c[2])
    assert np.max(np.abs(res)) < 0.001 * amp
    # 峰谷对称: E_max - E_min = 2A
    assert (pscan["E_kin_eV"].max() - pscan["E_kin_eV"].min()) ==         pytest.approx(2 * amp, rel=2e-3)


def test_pscan_crest_energy_vs_ref_particle(pscan):
    """PScan 相角 0 的能量 vs 同一运行 ref.001 末态参考粒子动能.

    该运行 Auto_Phase=F, 相角网格为绝对相位; ref 参考粒子携带
    ~2° 的注入相位偏移, 故允许 0.5% 偏差。
    """
    ref = read_ref_file(str(CAV / "astra_pscan"))
    e_ref = np.sqrt(ref.pz[-1] ** 2 + 0.51099895e6 ** 2) - 0.51099895e6
    i0 = int(np.argmin(np.abs(pscan["phase_deg"])))
    assert pscan["E_kin_eV"][i0] == pytest.approx(e_ref, rel=5e-3)


def test_pscan_compression_columns(pscan):
    """压缩因子/速度比列在 0..360° 全程有界且峰处 ~1."""
    assert np.all(np.isfinite(pscan["compression"]))
    assert np.all(pscan["compression"] > 0)
    i_crest = int(np.argmax(pscan["E_kin_eV"]))
    assert pscan["compression"][i_crest] == pytest.approx(1.0, abs=5e-3)
    assert pscan["beta_ratio"][i_crest] == pytest.approx(1.0, abs=5e-3)


# ---------------------------------------------------------------
# Scan 真实数据交叉验证
# ---------------------------------------------------------------
@pytest.fixture(scope="module")
def scan():
    return read_scan(DATA / "Example.Scan.001")


def test_scan_parameter_grid(scan):
    """para = 0.1..0.5 (S_min..S_max, S_numb=5), z = ZSTOP=1.5."""
    assert np.allclose(scan["para"], [0.1, 0.2, 0.3, 0.4, 0.5])
    assert np.allclose(scan["z"], 1.5)


def test_scan_fom_matches_golden_emit(scan):
    """FOM(1) 归一化水平发射度 [pi mm mrad] 与 golden Xemit 对照.

    golden 运行 MaxB=0.35 T, 落在扫描点 0.3 与 0.4 之间: 插值应
    再现 golden 末行 eps_nx; FOM(2)=束长 [mm], FOM(3)=平均能量
    [MeV] 与 Zemit 末行逐列相等。
    """
    xemit = np.loadtxt(DATA / "Example.Xemit.001")
    zemit = np.loadtxt(DATA / "Example.Zemit.001")
    eps_interp = np.interp(0.35, scan["para"], scan["FOM"][:, 0])
    assert eps_interp == pytest.approx(xemit[-1, 5], rel=5e-3)
    assert scan["FOM"][:, 1] == pytest.approx(zemit[-1, 3], rel=1e-4)
    assert scan["FOM"][:, 2] == pytest.approx(zemit[-1, 2], rel=1e-4)


# ---------------------------------------------------------------
# t 轴变体 / 速度 / 步长 (lineplot)
# ---------------------------------------------------------------
@pytest.fixture(scope="module")
def emit():
    return read_emit_files(str(DATA / "Example"), run="001")


def test_evolution_plots_time_axis(emit):
    import matplotlib.pyplot as plt
    from astra_tools.plot.emit_plots import (
        plot_envelope_evolution, plot_emittance_evolution,
        plot_bunch_length_evolution)
    for fn in (plot_envelope_evolution, plot_emittance_evolution,
               plot_bunch_length_evolution):
        fig = fn(emit, x_axis="t")
        assert "t [ns]" in fig.axes[0].get_xlabel()
        plt.close(fig)


def test_velocity_and_step_size_plots():
    import matplotlib.pyplot as plt
    from astra_tools.plot.emit_plots import (
        plot_velocity_evolution, plot_step_size_evolution)
    ref = read_ref_file(str(DATA / "Example"))
    fig = plot_velocity_evolution(ref)
    assert "velocity" in fig.axes[0].get_ylabel()
    plt.close(fig)
    fig = plot_step_size_evolution(ref)
    assert "step size [mm]" in fig.axes[0].get_ylabel()
    assert "z [m]" in fig.axes[0].get_xlabel()
    plt.close(fig)


# ---------------------------------------------------------------
# 孔径叠加 / 阴极表面场 / 激光 / 等离子体 (fieldplot)
# ---------------------------------------------------------------
def test_aperture_elements_from_namelist():
    from astra_tools.namelist.parse import parse_namelists
    from astra_tools.plot.advanced_plots import aperture_elements
    ap = parse_namelists(PROJECT_ROOT / "examples/Aperture/astra.in")["APERTURE"]
    els = aperture_elements(ap)
    assert len(els) == 2
    assert els[0]["r"] == pytest.approx(1.5e-3)      # Ap_R mm -> m
    assert els[0]["z1"] == pytest.approx(0.100)
    assert els[1]["r"] == pytest.approx(-50e-3)      # beam stop (negative)


def test_envelope_with_aperture_plot(emit):
    import matplotlib.pyplot as plt
    from astra_tools.plot.advanced_plots import (
        aperture_elements, plot_envelope_with_aperture)
    from astra_tools.namelist.parse import parse_namelists
    ap = parse_namelists(PROJECT_ROOT / "examples/Aperture/astra.in")["APERTURE"]
    fig = plot_envelope_with_aperture(emit, aperture_elements(ap))
    texts = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert any("aperture" in t for t in texts)
    assert "mm" in fig.axes[0].get_ylabel()
    plt.close(fig)


def test_cathode_emission_spch():
    import matplotlib.pyplot as plt
    from astra_tools.plot.advanced_plots import plot_cathode_emission
    fake = {"t": np.linspace(0, 1e-9, 5), "E_acc": np.ones(5) * 1e6,
            "E_spch": np.ones(5) * 1e5, "q": np.linspace(0, 1e-3, 5)}
    fig = plot_cathode_emission(fake, include_spch=True)
    texts = [t.get_text() for t in fig.axes[0].get_legend().get_texts()]
    assert any("spch" in t for t in texts)
    assert "t [ns]" in fig.axes[0].get_xlabel()
    plt.close(fig)


def test_laser_on_axis_plot():
    import matplotlib.pyplot as plt
    from astra_tools.plot.advanced_plots import plot_laser_on_axis
    fig = plot_laser_on_axis(PROJECT_ROOT / "examples/Cavity_Example/3D_test.ex",
                             unit="V/m")
    labels = fig.axes[0].get_xlabel() + fig.axes[1].get_xlabel()
    assert "z [mm]" in labels and "[ps]" in fig.axes[1].get_xlabel()
    plt.close(fig)


def test_plasma_profile_plot():
    import matplotlib.pyplot as plt
    from astra_tools.plot.advanced_plots import plot_plasma_profile
    fig = plot_plasma_profile(
        PROJECT_ROOT / "examples/Plasma_Example_1/PLASMA_flattop.txt",
        peak_density_cm3=1e17)
    assert "z [mm]" in fig.axes[0].get_xlabel()
    assert "arb. u." in fig.axes[0].get_ylabel()
    plt.close(fig)


# ---------------------------------------------------------------
# 核心电荷分数曲线 (postpro)
# ---------------------------------------------------------------
def test_central_charge_fraction_curves_full_bunch_matches_statistics():
    dist = read_distribution(DATA / "Example.0150.001")
    c = compute_central_charge_fraction_curves(dist, fractions=(1.0,))
    stats = compute_statistics(dist)
    assert c["sig_z"][0] == pytest.approx(stats.sig_z, rel=1e-9)
    assert c["emit_xn"][0] == pytest.approx(stats.emit_x_norm, rel=1e-9)
    assert c["n_particles"][0] == stats.n_active


def test_central_charge_fraction_plot():
    import matplotlib.pyplot as plt
    from astra_tools.plot.advanced_plots import plot_central_charge_fraction_curves
    dist = read_distribution(DATA / "Example.0150.001")
    fig = plot_central_charge_fraction_curves(dist)
    assert "core charge fraction" in fig.axes[0].get_xlabel()
    assert "mm mrad" in fig.axes[1].get_ylabel()
    plt.close(fig)

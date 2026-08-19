"""2026-08 对抗性审计修复: slice/发射度 (P1-5, F4, F3).

对应审计发现:
  P1  equi_energy 分箱是静默空操作 (退化为 z 等宽) -> 按能量等电荷分箱
  F4  正则动量螺线管修正不感知粒子种类 -> 正电荷种类符号翻转
  F3  Sigma eigen-发射度 x/y 标签按模长排序 -> 用本征向量识别平面
"""

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.constants import C_LIGHT, M_E_C2_EV
from astra_tools.distribution import Distribution
from astra_tools.analysis.statistics import compute_statistics
from astra_tools.analysis.slices import compute_slice_analysis
from astra_tools.io import read_sigma_file


def _chirped_beam(n=400, seed=5):
    """z 与能量强关联的束 (能啁啾): pz = 1e6 + 4e8*z."""
    rng = np.random.default_rng(seed)
    z = np.linspace(-1e-3, 1e-3, n) + rng.normal(0, 5e-6, n)
    pz = 1.0e6 + 4.0e8 * z
    return Distribution(
        x=rng.normal(0, 1e-4, n), y=rng.normal(0, 1e-4, n),
        z=z, px=np.zeros(n), py=np.zeros(n), pz=pz,
        clock=np.zeros(n), charge=np.full(n, -2e-3),
        status=np.full(n, 5),
        ref_momentum_eVc=1.0e6,
    )


def test_equi_energy_bins_by_energy_equal_charge():
    """equi_energy 必须按动能等电荷分箱, 与 equi_spaced 不同且每箱等 |q|."""
    dist = _chirped_beam()
    sa_e = compute_slice_analysis(dist, n_slices=10, binning="equi_energy")
    sa_z = compute_slice_analysis(dist, n_slices=10, binning="equi_spaced")
    # 与 z 等宽分箱输出不同 (旧实现逐位相同)
    assert not np.allclose(sa_e.z_edges, sa_z.z_edges)
    # z_edges 承载能量边界: 单调递增, 与束流动能范围一致
    from astra_tools.constants import kinetic_energy_from_momentum_vector
    e_all = kinetic_energy_from_momentum_vector(dist.px, dist.py, dist.pz)
    assert np.all(np.diff(sa_e.z_edges) > 0)
    assert sa_e.z_edges[0] == pytest.approx(float(np.min(e_all)), rel=1e-6)
    assert sa_e.z_edges[-1] == pytest.approx(float(np.max(e_all)), rel=1e-6)
    # 每箱 |q| 相等 (均匀宏电荷 -> 每箱粒子数在 n/n_slices 附近 ±2)
    counts = sa_e.n_particles[sa_e.n_particles > 0]
    assert np.ptp(counts) <= 2
    assert np.sum(counts) == dist.n_particle
    # 箱中心 (能量) 与各箱平均动能一致 (单调递增的啁啾)
    nz = sa_e.n_particles > 0
    assert np.all(np.diff(sa_e.mean_kinetic_energy_eV[nz]) > 0)
    assert np.all(sa_e.z_centers[nz] > 2e5)     # 中心是能量不是 z (z < 1 mm)


def test_equi_energy_plot_axis_is_energy():
    """equi_energy 的绘图横轴必须标能量, 不是 'z [mm]'."""
    import matplotlib.pyplot as plt
    from astra_tools.plot.slice_plots import plot_current_profile
    dist = _chirped_beam(n=300)
    sa = compute_slice_analysis(dist, n_slices=8, binning="equi_energy")
    fig = plot_current_profile(sa)
    try:
        assert "z [mm]" not in fig.axes[0].get_xlabel()
        assert "E" in fig.axes[0].get_xlabel()
    finally:
        plt.close("all")


def _solenoid_beam(species, n=2000, seed=9):
    """px = -c*Bz*y/2, py = +c*Bz*x/2 的束: 电子正则动量恒为零."""
    rng = np.random.default_rng(seed)
    x = rng.normal(0, 1e-3, n)
    y = rng.normal(0, 1e-3, n)
    bz = 0.5
    px = -0.5 * C_LIGHT * bz * y
    py = +0.5 * C_LIGHT * bz * x
    return Distribution(
        x=x, y=y, z=rng.normal(0, 1e-4, n),
        px=px, py=py, pz=np.full(n, 1.0e9),
        clock=np.zeros(n), charge=np.full(n, -2e-3),
        status=np.full(n, 5), index=np.full(n, species, dtype=np.int32),
        ref_momentum_eVc=1.0e9,
    )


def test_canonical_momentum_electron_sign_unchanged():
    """电子 (species=1): 正则动量消除 Bz 相关散角 -> sig_xp ≈ 0."""
    dist = _solenoid_beam(1)
    st = compute_statistics(dist, bz_on_axis_T=0.5)
    assert st.sig_xp < 0.05 * st.sig_yp or st.sig_xp < 1e-6
    assert st.sig_yp < 0.05 * st.sig_xp or st.sig_yp < 1e-6


def test_canonical_momentum_positron_sign():
    """正电子 (species=2): 修正符号必须翻转 -> sig_xp ≈ c*Bz*σy/p_ref."""
    dist = _solenoid_beam(2)
    st = compute_statistics(dist, bz_on_axis_T=0.5)
    expected = C_LIGHT * 0.5 * 1e-3 / 1.0e9   # c*Bz*σy/p_ref
    assert st.sig_xp == pytest.approx(expected, rel=0.1)


def test_sigma_eigen_emittance_plane_labels(tmp_path):
    """对角耦合块下, 大的 ε 属于 x 平面: 标签按本征向量而非模长."""
    mc = M_E_C2_EV
    sig11 = 1.0e-8          # σx² [m²]
    sig22 = 4.0e10 / mc**2  # px var -> εx = 20 m·eV/c
    sig33 = 1.0e-8          # σy²
    sig44 = 2.5e9 / mc**2   # py var -> εy = 5 m·eV/c (小!)
    sig55 = 1.0e-9
    sig66 = 2.5e11 / mc**2
    # 21 个上三角元 (row-major i<=j), 对角外全零
    upper = [sig11, 0, 0, 0, 0, 0,
             sig22, 0, 0, 0, 0,
             sig33, 0, 0, 0,
             sig44, 0, 0,
             sig55, 0,
             sig66]
    base = tmp_path / "test"
    with open(str(base) + ".Sigma.001", "w") as fh:
        fh.write("1.5 1000.0 " + " ".join("%.12g" % v for v in upper) + "\n")
    sd = read_sigma_file(str(base), "001")
    ex = 20.0 / mc
    ey = 5.0 / mc
    assert sd.enx[0] == pytest.approx(ex, rel=1e-6)
    assert sd.eny[0] == pytest.approx(ey, rel=1e-6)
    # 纵向: sqrt(sig55*sig66) 直接可验
    assert sd.enz[0] == pytest.approx(
        np.sqrt(1.0e-9 * 2.5e11), rel=1e-6)

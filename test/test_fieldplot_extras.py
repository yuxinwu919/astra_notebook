"""fieldplot 补充项测试 (F3 螺线管分量/四极 next page, F4 激光 rms 包络,
F5 电荷环位置/等离子体场 vs z/zeta).

对照手册 5.7 菜单 4/5 与 5.7.3 项 8/9。
"""

import numpy as np
import pytest

from astra_tools.plot.advanced_plots import (
    plot_curved_cathode_contour, plot_laser_envelope, plot_plasma_fields,
    plot_quadrupole_field)
from astra_tools.plot.field_plots import plot_solenoid_components

PROJECT_ROOT = __import__("pathlib").Path(__file__).resolve().parents[1]
CONTOUR = PROJECT_ROOT / "examples/Curved_Cathode_Example/Contour.dat"
PLASMA = PROJECT_ROOT / "examples/Plasma_Example_1/PLASMA_flattop.txt"


def _solenoid():
    from astra_tools.io.field_map import SolenoidField
    z = np.linspace(0, 0.3, 201)
    return SolenoidField(z=z, bz0=np.exp(-((z - 0.15) / 0.05) ** 2))


def _write_gaussian_laser(path, sigma0=3e-4, zc=1.5e-3, zr=0.8e-3):
    """写一个横向高斯、束腰在 zc 的合成 3D 激光场图 (紧凑头).

    File_A0 语义: 存 a₀⊥² = exp(-2r²/w²), w(z) = w0·sqrt(1+((z-zc)/zr)²),
    w0 = sigma0 为 1/e² 幅值半径。幅值 a₀ 分布 exp(-r²/w²) 的 rms = w/2。
    """
    nx = ny = 21
    nz = 31
    x = np.linspace(-1.5e-3, 1.5e-3, nx)
    y = np.linspace(-1.5e-3, 1.5e-3, ny)
    z = np.linspace(0, 3e-3, nz)
    XX, YY, ZZ = np.meshgrid(x, y, z, indexing="ij")
    w = sigma0 * np.sqrt(1.0 + ((ZZ - zc) / zr) ** 2)
    f = np.exp(-2.0 * (XX ** 2 + YY ** 2) / w ** 2)
    lines = ["%d %g %g" % (nx, x[0], x[1] - x[0]),
             "%d %g %g" % (ny, y[0], y[1] - y[0]),
             "%d %g %g" % (nz, z[0], z[1] - z[0])]
    lines.append(" ".join("%.8g" % v for v in f.ravel(order="F")))
    path.write_text("\n".join(lines))
    return z, sigma0, zc


def test_solenoid_components_renders():
    import matplotlib.pyplot as plt
    sol = _solenoid()
    fig = plot_solenoid_components(sol)
    try:
        titles = {ax.get_title() for ax in fig.axes}
        assert any("Bz alone" in t for t in titles)
        assert any("Br alone" in t for t in titles)
        assert any("radius" in t for t in titles)   # R3rd
        labels = [ax.get_ylabel() for ax in fig.axes]
        assert any("Bz" in s and "T" in s for s in labels)
        assert any("Br" in s and "T" in s for s in labels)
    finally:
        plt.close("all")


def test_quadrupole_two_cols_ideal():
    import matplotlib.pyplot as plt
    tmp = PROJECT_ROOT / "test" / "_quad_tmp.dat"
    try:
        z = np.linspace(0, 1.0, 50)
        np.savetxt(tmp, np.column_stack([z, np.full(50, 5.0)]))
        fig = plot_quadrupole_field(str(tmp))
        try:
            assert any("Gx" in ln.get_label() for ax in fig.axes
                       for ln in ax.lines)
            assert any("Gy" in ln.get_label() for ax in fig.axes
                       for ln in ax.lines)
        finally:
            plt.close("all")
    finally:
        tmp.unlink(missing_ok=True)


def test_quadrupole_four_cols_shows_bz():
    import matplotlib.pyplot as plt
    tmp = PROJECT_ROOT / "test" / "_quad_bz_tmp.dat"
    try:
        z = np.linspace(0, 1.0, 50)
        np.savetxt(tmp, np.column_stack([z, np.full(50, 5.0),
                                         np.full(50, -5.0),
                                         np.linspace(0, 0.1, 50)]))
        fig = plot_quadrupole_field(str(tmp))
        try:
            labels = [ln.get_label() for ax in fig.axes for ln in ax.lines]
            assert "Bz" in labels
        finally:
            plt.close("all")
    finally:
        tmp.unlink(missing_ok=True)


def test_laser_envelope_gaussian_focus(tmp_path):
    import matplotlib.pyplot as plt
    p = tmp_path / "laser.dat"
    z, sigma0, zc = _write_gaussian_laser(p)
    fig = plot_laser_envelope(str(p))
    try:
        # 对称: x/y 包络一致
        xl = [ln.get_ydata() for ln in fig.axes[0].lines
              if "x envelope" in ln.get_label()]
        yl = [ln.get_ydata() for ln in fig.axes[0].lines
              if "y envelope" in ln.get_label()]
        assert len(xl) == 1 and len(yl) == 1
        np.testing.assert_allclose(xl[0], yl[0], rtol=1e-3)
        # 包络在束腰 zc 处最小; a₀ 分布 exp(-r²/w²) 的 rms = w/√2,
        # 故 √(σx²+σy²) = w0 (1/e² 幅值半径)。旧 f² 权重得 w0/2,
        # 系统性偏小 (2026-08 审计 P1)。
        sig = np.sqrt(xl[0] ** 2 + yl[0] ** 2)
        kmin = int(np.argmin(sig))
        assert abs(z[kmin] - zc) < (z[1] - z[0]) * 2
        assert sig[kmin] == pytest.approx(sigma0 * 1e3, rel=0.05)
    finally:
        plt.close("all")


def test_laser_envelope_constant_gaussian(tmp_path):
    """恒定横向尺寸的激光场: rms 包络恒为常数."""
    import matplotlib.pyplot as plt
    p = tmp_path / "laser2.dat"
    z, sigma0, _ = _write_gaussian_laser(p, sigma0=2e-4, zr=1e9)  # 近似平行
    fig = plot_laser_envelope(str(p))
    try:
        xl = fig.axes[0].lines[0].get_ydata()
        # 去掉焦点附近 (zr 大, 几乎恒定), 全范围近常数
        assert np.ptp(xl) < 0.05 * np.mean(xl)
    finally:
        plt.close("all")


def test_cathode_rings_behind_surface(tmp_path):
    """电荷环位于阴极背面 (z 减小方向), 点数与表点数一致."""
    import matplotlib.pyplot as plt
    if not CONTOUR.exists():
        pytest.skip("Contour.dat 缺失")
    fig = plot_curved_cathode_contour(str(CONTOUR), show_rings=True)
    try:
        data = np.loadtxt(CONTOUR)
        z, tr = data[:, 0], data[:, 3]
        off = 0.1 * float(np.ptp(z))
        zr = z + off * (-tr)
        # 环在阴极背面: 比对应表面点更靠 -z (tr>0 恒成立)
        assert np.all(zr <= z + 1e-12)
        ring_line = [ln for ln in fig.axes[0].lines
                     if "rings" in ln.get_label()]
        assert len(ring_line) == 1
        assert len(ring_line[0].get_xdata()) == len(z)
    finally:
        plt.close("all")


def test_plasma_fields_z_and_zeta(tmp_path):
    import matplotlib.pyplot as plt
    if not PLASMA.exists():
        pytest.skip("PLASMA 文件缺失")
    for vs in ("z", "zeta"):
        fig = plot_plasma_fields(str(PLASMA), peak_density_cm3=1e17, vs=vs)
        try:
            ylabels = [ax.get_ylabel() for ax in fig.axes]
            assert any("Ez" in s for s in ylabels)
            assert any("density" in s for s in ylabels)
            if vs == "zeta":
                assert any("zeta" in ax.get_xlabel() for ax in fig.axes)
        finally:
            plt.close("all")


def test_plasma_fields_arb_units(tmp_path):
    """无峰值密度时 kp=1 (任意单位), 不抛异常且幅度归一."""
    import matplotlib.pyplot as plt
    if not PLASMA.exists():
        pytest.skip("PLASMA 文件缺失")
    fig = plot_plasma_fields(str(PLASMA), vs="z")
    try:
        ln = fig.axes[0].lines[0]
        assert np.all(np.abs(ln.get_ydata()) <= 1.0 + 1e-9)
    finally:
        plt.close("all")

"""2026-08 对抗性审计修复: 场展开与场图处理 (P1/P2/P3).

对应审计发现:
  P1   CavityField.field_at TM 展开缺 w²/c² 修正、多 r⁴ 项 (手册附录 I)
  P1   plot_laser_envelope 包络权重 f² -> 应 √f (laser.dat 存 a₀⊥²)
  P2   deck/solenoid.py 负 MaxB 静默跳过
  P3   螺线管仅 1 阶展开 (ASTRA 默认 3 阶) / z 非单调无校验 /
       3D 分量网格不一致直接拒绝 / E-B 族并存静默丢弃 / 就地改写无告警
"""

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.constants import C_LIGHT
from astra_tools.io.field_map import (
    CavityField, SolenoidField, read_cavity_field,
    read_3d_field_map_components, fix_laser_map_header,
)
from astra_tools.deck.solenoid import solenoid_bz_at_z

K = 2.0 * np.pi / 0.1          # k = 20π 1/m
OMEGA_SYNC = K * C_LIGHT       # 同步波: w = kc
W2C2_SYNC = (OMEGA_SYNC / C_LIGHT) ** 2


def _cos_field(z, k):
    return np.cos(k * z)


# ---------------------------------------------------------------- TM 展开

def test_tm_expansion_synchronous_cosine_manual():
    """同步纯余弦场: Ez 与 r 无关 (r² 项精确相消), 逐项对照手册附录 I."""
    z = np.linspace(-0.05, 0.05, 4001)
    cav = CavityField(z=z, ez0=_cos_field(z, K))
    r = 0.02
    zz = np.array([0.0, 0.1 / 8, -0.1 / 16, 0.0375])
    rr = np.full_like(zz, r)
    ez, er, bphi = cav.field_at(rr, zz, OMEGA_SYNC)
    np.testing.assert_allclose(ez, np.cos(K * zz), rtol=2e-5)
    np.testing.assert_allclose(er, (r / 2.0) * K * np.sin(K * zz), rtol=2e-5)
    np.testing.assert_allclose(
        bphi, (r / 2.0) * (OMEGA_SYNC / C_LIGHT**2) * np.cos(K * zz),
        rtol=2e-5)


def test_tm_expansion_nonsynchronous_cosine_manual():
    """非同步 (w != kc): w²/c² 修正项不可忽略, 逐项对照手册附录 I."""
    omega = 1.5 * K * C_LIGHT
    w2c2 = (omega / C_LIGHT) ** 2
    z = np.linspace(-0.05, 0.05, 4001)
    cav = CavityField(z=z, ez0=_cos_field(z, K))
    r = 0.02
    zz = np.array([0.0, 0.1 / 8, -0.1 / 16])
    rr = np.full_like(zz, r)
    ez, er, bphi = cav.field_at(rr, zz, omega)
    c0 = np.cos(K * zz)
    s0 = np.sin(K * zz)
    ez_exp = c0 - (r**2 / 4.0) * (w2c2 - K**2) * c0
    er_exp = (-(r / 2.0) * (-K * s0)
              + (r**3 / 16.0) * (K**3 * s0 + w2c2 * (-K * s0)))
    bp_exp = ((r / 2.0) * c0
              - (r**3 / 16.0) * (w2c2 - K**2) * c0) * (omega / C_LIGHT**2)
    np.testing.assert_allclose(ez, ez_exp, rtol=2e-5)
    np.testing.assert_allclose(er, er_exp, rtol=2e-5)
    np.testing.assert_allclose(bphi, bp_exp, rtol=2e-5)


def test_tm_expansion_static_cosine_no_r4_term():
    """静场 (w=0) 的 TM 展开: 手册无 r⁴ 项 (旧代码有, 超 3 阶)."""
    z = np.linspace(-0.05, 0.05, 4001)
    cav = CavityField(z=z, ez0=_cos_field(z, K))
    r = 0.03
    zz = np.array([0.0, 0.1 / 8])
    rr = np.full_like(zz, r)
    ez, er, bphi = cav.field_at(rr, zz, 0.0)
    c0 = np.cos(K * zz)
    ez_exp = c0 + (r**2 / 4.0) * K**2 * c0        # 无 r⁴ 项
    np.testing.assert_allclose(ez, ez_exp, rtol=2e-5)
    assert np.all(bphi == 0.0)


# ---------------------------------------------------------------- 螺线管

def test_solenoid_third_order_manual():
    """螺线管 3 阶展开 (ASTRA S_higher_order 默认): 含 r⁴/64 与 r³/16 项."""
    z = np.linspace(-0.05, 0.05, 4001)
    sol = SolenoidField(z=z, bz0=_cos_field(z, K))
    r = 0.02
    zz = np.array([0.0, 0.1 / 8])
    rr = np.full_like(zz, r)
    br, bz = sol.field_at(rr, zz)
    c0 = np.cos(K * zz)
    s0 = np.sin(K * zz)
    bz_exp = c0 + (r**2 / 4.0) * K**2 * c0 + (r**4 / 64.0) * K**4 * c0
    br_exp = (r / 2.0) * K * s0 + (r**3 / 16.0) * K**3 * s0
    np.testing.assert_allclose(bz, bz_exp, rtol=2e-5)
    np.testing.assert_allclose(br, br_exp, rtol=2e-5)


# ---------------------------------------------------------------- 读取校验

def test_read_cavity_field_rejects_nonmonotonic_z(tmp_path):
    """z 非单调的场表必须报错 (np.interp/gradient 会静默给错)."""
    p = tmp_path / "cav.dat"
    p.write_text("0.0 1.0\n0.1 2.0\n0.05 3.0\n0.2 4.0\n")
    with pytest.raises(ValueError):
        read_cavity_field(p)


# ---------------------------------------------------------------- 3D 分量

def _write_compact_3d(path, nx, ny, nz, z0=0.0, dz=1e-3, val=1.0):
    x = np.linspace(-1e-3, 1e-3, nx)
    y = np.linspace(-1e-3, 1e-3, ny)
    z = np.linspace(z0, z0 + (nz - 1) * dz, nz)
    f = np.full((nx, ny, nz), val)
    lines = ["%d %g %g" % (nx, x[0], x[1] - x[0]),
             "%d %g %g" % (ny, y[0], y[1] - y[0]),
             "%d %g %g" % (nz, z[0], dz)]
    lines.append(" ".join("%.8g" % v for v in f.ravel(order="F")))
    path.write_text("\n".join(lines))


def test_3d_components_different_grids_interpolate(tmp_path):
    """各分量网格不同 (手册允许): 插值到首分量网格 + 告警, 不拒绝."""
    base = tmp_path / "map"
    _write_compact_3d(Path(str(base) + ".ex"), 3, 3, 4, val=1.0)
    _write_compact_3d(Path(str(base) + ".ey"), 3, 3, 5, val=2.0)
    _write_compact_3d(Path(str(base) + ".ez"), 3, 3, 4, val=3.0)
    with pytest.warns(UserWarning):
        m = read_3d_field_map_components(base)
    assert m.fx.shape == (3, 3, 4)
    assert m.fz.shape == (3, 3, 4)          # 同网格分量原样
    assert m.fy.shape == (3, 3, 4)          # 已插值到公共网格
    np.testing.assert_allclose(m.fx, 1.0)
    np.testing.assert_allclose(m.fy, 2.0, rtol=1e-9)


def test_3d_components_eb_families_warn(tmp_path):
    """E/B 两族并存: 只显示 E 族但必须告警 (不再静默丢弃 B 族)."""
    base = tmp_path / "map2"
    _write_compact_3d(Path(str(base) + ".ex"), 3, 3, 4, val=1.0)
    _write_compact_3d(Path(str(base) + ".bx"), 3, 3, 4, val=0.5)
    with pytest.warns(UserWarning):
        m = read_3d_field_map_components(base)
    assert m.quantity == "E"
    assert m.unit == "V/m"


# ---------------------------------------------------------------- 负 MaxB

def test_negative_maxB_solenoid_not_skipped(tmp_path):
    """负 MaxB (反向螺线管) 是合法输入: 按符号缩放并告警, 不静默跳过."""
    d = tmp_path
    (d / "f.dat").write_text(
        "\n".join("%g %g" % (z, 1.0) for z in np.linspace(0, 2, 101)))
    (d / "astra.in").write_text(
        "&SOLENOID\n LBField=T,\n File_Bfield(1)='f.dat',\n"
        " MaxB(1)=-0.2,\n S_pos(1)=0.0,\n/\n")
    with pytest.warns(UserWarning):
        bz = solenoid_bz_at_z(d / "astra.in", 1.0)
    assert bz is not None
    assert bz == pytest.approx(-0.2)         # 符号保留 (反向场)


# ---------------------------------------------------------------- laser

def test_fix_laser_map_header_warns_inplace(tmp_path):
    """就地改写用户文件必须告警."""
    p = tmp_path / "laser.dat"
    p.write_text("8.1e+01 0.0 1e-3\n9.0e+00 0.0 2e-3\n4.0e+01 0.0 3e-3\n0.1\n")
    with pytest.warns(UserWarning):
        fix_laser_map_header(p)


def test_laser_envelope_weight_is_amplitude(tmp_path):
    """laser.dat 存 a₀⊥²: 包络权重必须 a₀=√f; f² 权重使 rms 偏小 1/√2."""
    from astra_tools.plot.advanced_plots import plot_laser_envelope
    import matplotlib.pyplot as plt
    p = tmp_path / "laser.dat"
    # a₀² = exp(-2 r² / w²), w = 1/e² 幅值半径, 束腰 w0=3e-4 (平行束)
    nx = ny = 41
    nz = 31
    w0 = 3e-4
    x = np.linspace(-1.5e-3, 1.5e-3, nx)
    y = np.linspace(-1.5e-3, 1.5e-3, ny)
    z = np.linspace(0, 3e-3, nz)
    XX, YY = np.meshgrid(x, y, indexing="ij")
    f = np.exp(-2.0 * (XX**2 + YY**2) / w0**2)[:, :, None] * np.ones((1, 1, nz))
    lines = ["%d %g %g" % (nx, x[0], x[1] - x[0]),
             "%d %g %g" % (ny, y[0], y[1] - y[0]),
             "%d %g %g" % (nz, z[0], z[1] - z[0])]
    lines.append(" ".join("%.8g" % v for v in f.ravel(order="F")))
    p.write_text("\n".join(lines))
    fig = plot_laser_envelope(str(p))
    try:
        xl = [ln.get_ydata() for ln in fig.axes[0].lines
              if "x envelope" in ln.get_label()][0]
        # a₀ 分布 exp(-r²/w²) 的 rms = w/√2 (f² 权重则得 w/(2√2), 偏小 √2)
        assert np.mean(xl[10:-10]) * 1e-3 == pytest.approx(w0 / np.sqrt(2.0), rel=1e-3)
    finally:
        plt.close("all")

"""TE 模场与场展开半径 R3rd 的单元测试 (手册 6.9 / 8 章).

R3rd 公式 (手册 8 章):
  TM (Ez0):  R3rd_Er = sqrt(|0.08 Ez'| / |Ez''' + w²/c² Ez'|)
             R3rd_Bφ = sqrt(|0.08 Ez0| / |Ez'' + w²/c² Ez0|)
  TE (Bz0):  R3rd_Eφ = sqrt(|0.08 Bz0| / |Bz'' + w²/c² Bz0|)
             R3rd_Br = sqrt(|0.08 Bz'| / |Bz'' + w²/c² Bz'|)
  螺线管:    R3rd = sqrt(|0.08 Bz'| / |Bz'''|)

TE 离轴展开 (手册 8 章):
  Bz(r) = Bz0 - (r²/4)(Bz0'' + w²/c² Bz0)
  Br(r) = -(r/2) Bz0' + (r³/16)(Bz0''' + w²/c² Bz0')
  Eφ(r) = [(r/2) Bz0 - (r³/16)(Bz0'' + w²/c² Bz0)] ω
"""

import numpy as np
import pytest

from astra_tools.io.field_map import CavityField, SolenoidField, TEField

C_LIGHT = 299792458.0
OMEGA = 2 * np.pi * 1.3e9            # 1.3 GHz L-band
W2C2 = (OMEGA / C_LIGHT) ** 2        # ~741.3 m^-2


def _gaussian(n=301, zmin=0.0, zmax=0.3, mu=0.15, sigma=0.05, amp=1.0):
    z = np.linspace(zmin, zmax, n)
    return z, amp * np.exp(-((z - mu) ** 2) / (2 * sigma ** 2))


def test_te_field_at_constant_bz0_manual():
    """Bz0 恒定: Br=0, Bz=B0(1-r²w²/4c²), Eφ=B0 ω(r/2 - r³w²/16c²)."""
    z = np.linspace(0, 0.3, 51)
    te = TEField(z=z, bz0=np.full_like(z, 2.0))
    r = np.array([0.0, 0.005, 0.01, 0.02])
    zz = np.array([0.15])
    bz, br, ephi = te.field_at(r, zz, OMEGA)
    expected_bz = 2.0 * (1.0 - r ** 2 / 4.0 * W2C2)
    expected_br = np.zeros_like(r)
    expected_ephi = 2.0 * OMEGA * (r / 2.0 - r ** 3 / 16.0 * W2C2)
    np.testing.assert_allclose(bz, expected_bz, rtol=1e-6)
    np.testing.assert_allclose(br, expected_br, atol=1e-12)
    np.testing.assert_allclose(ephi, expected_ephi, rtol=1e-6)
    # 轴上: Bz=B0, Br=Eφ=0
    bz0a, br0a, ep0a = te.field_at(np.array([0.0]), zz, OMEGA)
    assert bz0a[0] == pytest.approx(2.0)
    assert br0a[0] == pytest.approx(0.0, abs=1e-12)
    assert ep0a[0] == pytest.approx(0.0, abs=1e-12)


def test_te_field_at_linear_ramp_manual():
    """Bz0(z)=B0+cz (线性): 各阶导解析, 核对数值展开."""
    z = np.linspace(0, 0.3, 501)
    b0, c = 1.0, 3.0
    te = TEField(z=z, bz0=b0 + c * z)
    r = np.array([0.01, 0.02])
    zz = np.array([0.1, 0.2])
    bz, br, ephi = te.field_at(r, zz, OMEGA)
    b0a, cp = b0 + c * zz, c
    ebz = b0a - (r ** 2 / 4.0) * (0.0 + W2C2 * b0a)
    ebr = -(r / 2.0) * cp + (r ** 3 / 16.0) * (0.0 + W2C2 * cp)
    eep = ((r / 2.0) * b0a - (r ** 3 / 16.0) * (0.0 + W2C2 * b0a)) * OMEGA
    np.testing.assert_allclose(bz, ebz, rtol=1e-5)
    np.testing.assert_allclose(br, ebr, rtol=1e-5)
    np.testing.assert_allclose(ephi, eep, rtol=1e-5)


def test_te_expansion_radius_constant_bz0():
    """恒定 Bz0: R3rd_Eφ = c/ω·sqrt(0.08), R3rd_Br 退化 (~0, 噪声级)."""
    z = np.linspace(0, 0.3, 101)
    te = TEField(z=z, bz0=np.full_like(z, 1.5))
    r_ep, r_br = te.expansion_radius(OMEGA)
    expected = np.sqrt(0.08) * C_LIGHT / OMEGA
    np.testing.assert_allclose(r_ep[1:-1], expected, rtol=1e-5)
    # 无 Bz 导数 -> Br 半径退化: 应 NaN 或梯度浮点噪声级小值。
    # 退化分量无物理意义, 只要求不产生异常大的有限值。
    finite = r_br[np.isfinite(r_br)]
    assert finite.size == 0 or float(finite.max()) < 2 * expected


def test_cavity_expansion_radius_constant_ez0():
    """恒定 Ez0: R3rd_Bφ = c/ω·sqrt(0.08), R3rd_Er 退化 (噪声级小值)."""
    z = np.linspace(0, 0.3, 101)
    cav = CavityField(z=z, ez0=np.full_like(z, 10.0))
    r_er, r_bp = cav.expansion_radius(OMEGA)
    expected = np.sqrt(0.08) * C_LIGHT / OMEGA
    np.testing.assert_allclose(r_bp[1:-1], expected, rtol=1e-5)
    # 无 Ez 导数 -> Er 半径退化 (梯度浮点噪声放大, 应远小于主半径)
    finite = r_er[np.isfinite(r_er)]
    assert finite.size == 0 or float(finite.max()) < expected


def test_cavity_expansion_radius_manual_formula():
    """TM: R3rd_Er = sqrt(|0.08 Ez'|/|Ez'''+w²c²Ez'|) 独立重算核对."""
    z, ez0 = _gaussian()
    cav = CavityField(z=z, ez0=ez0)
    r_er, r_bp = cav.expansion_radius(OMEGA)
    dez = np.gradient(ez0, z)
    d2ez = np.gradient(dez, z)
    d3ez = np.gradient(d2ez, z)
    with np.errstate(divide="ignore", invalid="ignore"):
        exp_er = np.sqrt(np.abs(0.08 * dez) / np.abs(d3ez + W2C2 * dez))
        exp_bp = np.sqrt(np.abs(0.08 * ez0) / np.abs(d2ez + W2C2 * ez0))
    np.testing.assert_allclose(r_er, exp_er, rtol=1e-9)
    np.testing.assert_allclose(r_bp, exp_bp, rtol=1e-9)
    # 区域内应基本为正且有限 (除退化为 0 的点)
    ok = np.isfinite(r_er) & np.isfinite(r_bp)
    assert ok.sum() > len(z) // 2


def test_solenoid_expansion_radius_manual_formula():
    """螺线管: R3rd = sqrt(|0.08 Bz'|/|Bz'''|) 独立重算核对."""
    z, bz0 = _gaussian()
    sol = SolenoidField(z=z, bz0=bz0)
    r = sol.expansion_radius()
    dbz = np.gradient(bz0, z)
    d3bz = np.gradient(np.gradient(dbz, z), z)
    with np.errstate(divide="ignore", invalid="ignore"):
        exp = np.sqrt(np.abs(0.08 * dbz) / np.abs(d3bz))
    np.testing.assert_allclose(r, exp, rtol=1e-9)
    ok = np.isfinite(r)
    assert ok.sum() > len(z) // 2


def test_read_te_field(tmp_path):
    """read_te_field: 两列表 (z, Bz) 解析, 数据往返一致."""
    from astra_tools.io import read_te_field
    p = tmp_path / "TE_test.dat"
    z = np.linspace(0, 0.3, 21)
    bz0 = np.sin(z * 20)
    np.savetxt(p, np.column_stack([z, bz0]))
    te = read_te_field(p)
    assert isinstance(te, TEField)
    assert te.peak_field_arb == pytest.approx(float(np.max(np.abs(bz0))))
    np.testing.assert_allclose(te.z, z)
    np.testing.assert_allclose(te.bz0, bz0)
    # scaled: 峰值归一
    te2 = te.scaled(0.5)
    assert te2.peak_field_arb == pytest.approx(0.5)
    np.testing.assert_allclose(te2.bz0, bz0 * (0.5 / te.peak_field_arb))

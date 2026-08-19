"""2026-08 对抗性审计修复: BFF 中性束回退与时间坐标测试增强.

对应审计发现:
  P3-1 近中性束团 BFF 静默归零 -> |q| 归一化结构因子回退 + 告警
  P2   时间坐标测试是恒等式 -> 真实数据 σ_t = σ_z/(βc) 对照
  P2   BFF 无独立解析对照 -> 双点电荷 cos²(kd/2) 解析对照
"""

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.analysis.bff import compute_bff
from astra_tools.analysis.time import bunch_time
from astra_tools.io import read_distribution
from astra_tools.distribution import Distribution
from astra_tools.constants import C_LIGHT, beta_from_gamma, gamma_from_momentum

DATA = PROJECT_ROOT / "examples" / "Manual_Example"


def test_neutral_bunch_bff_q_abs_fallback():
    """近中性束 (Σq=0): |q| 归一化结构因子回退 + 告警, 不再静默归零.

    两点电荷 ±q 相距 d: F̃(k) = cos(kd/2) -> BFF = cos²(kd/2),
    BFF(0) = 1。旧实现返回全零, 丢失 k≠0 处的结构因子信息。
    """
    d = 1e-3
    z = np.array([-d / 2.0, d / 2.0])
    q = np.array([1.0, -1.0])
    k = np.array([0.0, 1e3, 2e3, 3e3])
    with pytest.warns(UserWarning):
        b = compute_bff(z, q, kmin=0.0, kmax=3e3, nk=4,
                        log_spaced=False, method="direct")
    expected = np.cos(k * d / 2.0) ** 2
    assert b.bff[0] == pytest.approx(1.0, rel=1e-12)
    np.testing.assert_allclose(b.bff, expected, rtol=1e-12, atol=1e-12)


def test_bff_two_charges_analytic_direct_and_fft():
    """双等电荷: F(k) = cos(kd/2), BFF = cos²(kd/2) — 独立解析对照."""
    d = 2e-3
    z = np.array([-d / 2.0, d / 2.0])
    q = np.array([1.0, 1.0])
    for method in ("direct", "fft"):
        b = compute_bff(z, q, kmin=10.0, kmax=2e3, nk=50,
                        log_spaced=True, method=method)
        expected = np.cos(b.k * d / 2.0) ** 2
        np.testing.assert_allclose(b.bff, expected, rtol=2e-3, atol=1e-3)


def test_bunch_time_real_data_matches_z_over_beta_c():
    """真实 ASTRA 输出: active 粒子时间坐标 σ_t = σ_z/(β̄c) (1% 容差)."""
    dist = read_distribution(DATA / "Example.0150.001")
    t = bunch_time(dist)
    m = dist.active
    z_avr = float(np.mean(dist.z[m]))
    pz_avr = float(np.mean(dist.pz[m]))
    beta_bar = beta_from_gamma(gamma_from_momentum(pz_avr))
    sigma_z = float(np.std(dist.z[m] - z_avr))
    sigma_t = float(np.std(t[m]))
    assert sigma_t == pytest.approx(sigma_z / (beta_bar * C_LIGHT), rel=0.01)


def test_bunch_time_not_started_uses_clock_branch():
    """未发射粒子 (status -1..-6) 的时间坐标 = clock (不告警时)."""
    n = 10
    clock = np.linspace(1e-12, 10e-12, n)
    dist = Distribution(
        x=np.zeros(n), y=np.zeros(n), z=np.zeros(n),
        px=np.zeros(n), py=np.zeros(n), pz=np.full(n, 1.0e6),
        clock=clock, charge=np.full(n, -2e-3),
        status=np.full(n, -1), ref_momentum_eVc=1.0e6,
    )
    t = bunch_time(dist, warn_mixed=False)
    np.testing.assert_allclose(t, clock, rtol=1e-12)

"""2026-08 对抗性审计修复: 动能必须用全动量 |p| 而非仅 pz (P2-2).

大发散束 (如 90deg 探针束) 下 pz-only 的动能被低估可达 60%;
修复后 statistics / slices / cuts 的每粒子动能统一为
E_kin = sqrt((px^2+py^2+pz^2) + m^2c^4) - mc^2。
"""

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.constants import M_E_C2_EV
from astra_tools.distribution import Distribution
from astra_tools.analysis.statistics import compute_statistics
from astra_tools.analysis.slices import compute_slice_analysis
from astra_tools.analysis.cuts import cut_distribution


def _divergent_beam(n=200, seed=3):
    """pz=1e6 eV/c、px=3e6 eV/c、py=0 的大发散束 (90deg 场景)."""
    rng = np.random.default_rng(seed)
    return Distribution(
        x=rng.normal(0, 1e-4, n), y=rng.normal(0, 1e-4, n),
        z=np.linspace(-1e-3, 1e-3, n),
        px=np.full(n, 3.0e6), py=np.zeros(n), pz=np.full(n, 1.0e6),
        clock=np.zeros(n), charge=np.full(n, -2e-3),
        status=np.full(n, 5),
        ref_momentum_eVc=1.0e6,
    )


def _e_full(px, py, pz):
    return np.sqrt(px**2 + py**2 + pz**2 + M_E_C2_EV**2) - M_E_C2_EV


def test_statistics_mean_energy_uses_full_momentum():
    """mean_E_kin 必须 = E_kin(|p|), 而非 E_kin(pz)."""
    dist = _divergent_beam()
    stats = compute_statistics(dist)
    expected = float(_e_full(3.0e6, 0.0, 1.0e6))
    assert stats.mean_E_kin_eV == pytest.approx(expected, rel=1e-9)
    # 全粒子动能相同 -> 能散为零
    assert stats.sig_E_eV == pytest.approx(0.0, abs=1e-6)
    assert stats.sig_E_over_E == pytest.approx(0.0, abs=1e-12)


def test_slice_energy_and_beta_use_full_momentum():
    """slice 平均动能与纵向速度 beta 都用全动量."""
    dist = _divergent_beam(n=300)
    sa = compute_slice_analysis(dist, n_slices=5)
    expected = float(_e_full(3.0e6, 0.0, 1.0e6))
    nz = sa.n_particles > 0
    assert np.all(sa.mean_kinetic_energy_eV[nz] == pytest.approx(expected, rel=1e-9))
    # v_z/c = pz/sqrt(p^2 + m^2c^4) < beta(gamma(mean pz)) 对大发散束
    e_tot = np.sqrt(3.0e6**2 + 1.0e6**2 + M_E_C2_EV**2)
    vz_c = 1.0e6 / e_tot
    assert np.all(sa.beta_per_slice[nz] == pytest.approx(vz_c, rel=1e-9))


def test_cut_energy_window_uses_full_momentum():
    """e_range 按全动量动能切割: pz-only 值 (0.61 MeV) 的窗口应切不到粒子."""
    dist = _divergent_beam(n=100)
    _, mask_none = cut_distribution(dist, e_range=(2.6e6, 2.8e6))   # 含真动能
    assert not np.any(mask_none)                                    # 全保留
    _, mask_all = cut_distribution(dist, e_range=(0.5e6, 0.7e6))    # pz-only 值
    assert np.all(mask_all)                                         # 全切

"""R2a 低能/大发散交叉验证: 低能束 (5 MeV 与 gamma≈1.01) 真跑 golden。

第三阶段 Task 2: 既有交叉验证 golden 全部为 ~1 GeV 低发散束; 本测试
补齐低能区 (γ≈1.01 ~ γ≈10.8) 的 γβ 归一化口径验证。金样由本地 ASTRA
V4.0 (macOS Apple Silicon) 真跑产生 (deck: examples/LowEnergy_Validation/
LowEnergy.in, 纯漂移、无螺线管、无空间电荷, ZSTOP=0.5 m):

  * LowEnergy.Xemit/Zemit.001  — 5 MeV 束 (γ=10.78, 5.487 MeV/c)
  * LowEnergy.0050.001         — 同 run 相位 dump (z=0.5 m)
  * LowEnergy_lowg.*          — E_kin=5.1 keV 束 (γ=1.00998, 72.4 keV/c)

对照口径与 test_cross_validation.py 一致 (容差同源): Xemit/Zemit 末行
vs compute_statistics(相位 dump)。重点断言: γ≈1.01 下 γβ 归一化发射度
与 ASTRA 一致 — 验证 "γ 从动量算" (gamma = sqrt(1+(p/mc)^2)) 口径在
低能区成立; 判别力断言: 若把 p_ref 误当动能 (gamma = 1+p/mc^2),
βγ 会大 ~3.9 倍, 测试必须能区分 (防退化)。
"""

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.io import read_distribution
from astra_tools.io.astra_emit import parse_output_file
from astra_tools.analysis.statistics import compute_statistics
from astra_tools.constants import M_E_C2_EV, gamma_from_momentum

DATA = PROJECT_ROOT / "examples" / "LowEnergy_Validation"
GOLDEN = DATA / "golden"

# 与 test_cross_validation.py 同源容差
REL_MEAN_X = 5e-3
REL_SIG_X = 5e-4
REL_SIG_XP = 5e-4
REL_EPS_NX = 5e-3
REL_MEAN_E = 1e-4
REL_SIG_Z = 1e-3
REL_SIG_E = 1e-3
REL_EPS_ZN = 2e-3


def _xemit_stats(name: str):
    """read phase dump + compute_statistics, 返回 (stats, xemit_last, zemit_last)."""
    dist = read_distribution(GOLDEN / (name + ".0050.001"))
    stats = compute_statistics(dist)   # 无螺线管 -> bz=0
    xemit = parse_output_file(GOLDEN / (name + ".Xemit.001"))
    zemit = parse_output_file(GOLDEN / (name + ".Zemit.001"))
    return stats, xemit, zemit


def _check_transverse(stats, xemit):
    """Xemit 末行 vs compute_statistics (口径同 test_cross_validation.py)."""
    last = xemit["mean_z"][-1]
    assert last == pytest.approx(0.5, abs=1e-4)     # 样本确在 z=0.5 m
    assert stats.mean_x * 1e3 == pytest.approx(xemit["mean_x"][-1] * 1e3,
                                               rel=REL_MEAN_X)
    assert stats.sig_x * 1e3 == pytest.approx(xemit["sigma_x"][-1] * 1e3,
                                              rel=REL_SIG_X)
    assert stats.sig_xp * 1e3 == pytest.approx(xemit["sigma_xp"][-1] * 1e3,
                                               rel=REL_SIG_XP)
    # Xemit 第 6 列 = eps_n 单位 1e-6 m.rad (mm mrad)
    assert stats.emit_x_norm * 1e6 == pytest.approx(
        xemit["norm_emit_x"][-1] * 1e6, rel=REL_EPS_NX)


def _check_longitudinal(stats, zemit):
    """Zemit 末行 vs compute_statistics."""
    last = zemit["mean_z"][-1]
    assert last == pytest.approx(0.5, abs=1e-4)
    assert stats.mean_E_kin_eV == pytest.approx(zemit["mean_kinetic_energy"][-1],
                                                rel=REL_MEAN_E)
    assert stats.sig_z * 1e3 == pytest.approx(zemit["sigma_z"][-1] * 1e3,
                                              rel=REL_SIG_Z)
    assert stats.sig_E_eV == pytest.approx(zemit["sigma_energy"][-1],
                                           rel=REL_SIG_E)
    # 纵向 RMS 发射度 (含 <zE> 协方差): ASTRA 以 eV.m (= keV.mm 数值) 存
    assert stats.emit_z_eVm == pytest.approx(zemit["norm_emit_z"][-1],
                                             rel=REL_EPS_ZN)


def test_5mev_transverse_matches_xemit():
    stats, xemit, _ = _xemit_stats("LowEnergy")
    _check_transverse(stats, xemit)


def test_5mev_longitudinal_matches_zemit():
    stats, _, zemit = _xemit_stats("LowEnergy")
    _check_longitudinal(stats, zemit)


def test_5mev_reference_momentum():
    dist = read_distribution(GOLDEN / "LowEnergy.0050.001")
    # 参考动量 ~5.487 MeV/c <-> E_kin 5.0 MeV
    assert dist.ref_momentum_eVc == pytest.approx(5.4873e6, rel=1e-4)
    g = gamma_from_momentum(dist.ref_momentum_eVc)
    assert g == pytest.approx(10.7848, rel=1e-4)


def test_lowgamma_gamma1p01_normalized_emittance_matches_astra():
    """核心断言: γ≈1.01 下 γβ 归一化发射度与 ASTRA 一致.

    E_kin=5.1 keV -> p=72.4 keV/c, γ=1.00998, βγ=0.1416。ASTRA Xemit
    eps_n = 0.1362 mm mrad; 若按错误口径 γ = 1 + p_ref/mc^2 (把动量当
    动能) 则 βγ=0.551, eps_n 会大 3.9 倍 — 判别力断言保证本测试对
    "γ 从动量算" 口径失效立即失败。
    """
    stats, xemit, _ = _xemit_stats("LowEnergy_lowg")
    _check_transverse(stats, xemit)

    # γ 从动量算 (本项目口径)
    assert stats.gamma == pytest.approx(1.009980, rel=1e-4)
    assert stats.beta_rel == pytest.approx(0.140236, rel=1e-3)
    bg = np.sqrt(max(stats.gamma**2 - 1.0, 0.0))
    assert bg == pytest.approx(0.141635, rel=1e-4)

    # 判别力: 把 p_ref 误当动能 -> βγ 大 ~3.9 倍, 与 ASTRA 显著不符
    g_wrong = 1.0 + stats.ref_momentum_eVc / M_E_C2_EV
    bg_wrong = np.sqrt(max(g_wrong**2 - 1.0, 0.0))
    assert bg_wrong / bg > 3.0
    eps_wrong = bg_wrong * stats.emit_x_geom * 1e6
    assert abs(eps_wrong - xemit["norm_emit_x"][-1] * 1e6) \
        / (xemit["norm_emit_x"][-1] * 1e6) > 2.0


def test_lowgamma_longitudinal_matches_zemit():
    stats, _, zemit = _xemit_stats("LowEnergy_lowg")
    _check_longitudinal(stats, zemit)
    assert stats.mean_E_kin_eV == pytest.approx(5104.4, rel=1e-3)


def test_lowgamma_reference_momentum():
    dist = read_distribution(GOLDEN / "LowEnergy_lowg.0050.001")
    assert dist.ref_momentum_eVc == pytest.approx(7.2375e4, rel=1e-4)
    g = gamma_from_momentum(dist.ref_momentum_eVc)
    assert g == pytest.approx(1.00998, rel=1e-4)

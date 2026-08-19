"""R2b 90deg 大发散探针束交叉验证: ASTRA Zemit sigma_E 口径判定.

科学结论 (本项目笔记 07 的关键事实): **ASTRA Zemit 的每粒子动能
E_kin 与 sigma_E 采用全动量 |p| 口径**
(E_kin = sqrt(px^2+py^2+pz^2 + m^2c^4) - mc^2), 与本项目
compute_statistics 的 2026-08 全动量修复一致。

证据: examples/90deg_bend_Example/golden/Section1.Zemit.001 真跑输出
(deck Section1.in + EmitS=T; 与归档相位 dump Section1.0100.001 同 run,
已验字节一致)。探针束经 90deg 二极偏转后 px≈3.5 MeV/c、pz≈1.87 MeV/c,
|p|=3.975 MeV/c — 大发散判别力极强:

  口径            mean_E [MeV]  sig_E [keV]   eps_z [eV.m]
  ASTRA Zemit      3.4967       21.180         94.872
  全动量 |p|       3.496713     21.1800        94.872    <- 吻合 (< 1e-5)
  pz-only          1.425634     27.7626       126.287    <- 差 59%/31%/33%

参考粒子 pz=1.868 MeV/c -> Xemit 的 eps_n 用 gamma(p_ref) 归一化
(bg=3.656), 同样与本实现一致。
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
from astra_tools.constants import M_E_C2_EV

DATA = PROJECT_ROOT / "examples" / "90deg_bend_Example"
GOLDEN = DATA / "golden"
PHASE = GOLDEN / "Section1.0100.001"
ZEMIT = GOLDEN / "Section1.Zemit.001"
XEMIT = GOLDEN / "Section1.Xemit.001"


@pytest.fixture(scope="module")
def probe():
    """90deg 探针束相位 dump (归档 golden, 与 Zemit 同 run)."""
    return read_distribution(PHASE)


def _e_full(px, py, pz):
    return np.sqrt(px**2 + py**2 + pz**2 + M_E_C2_EV**2) - M_E_C2_EV


def _e_pz(pz):
    return np.sqrt(pz**2 + M_E_C2_EV**2) - M_E_C2_EV


def test_probe_beam_is_large_divergence():
    """探针束确为大发散: px/pz ~ 1.9, |p| 与 pz 差异巨大 (判别力前提)."""
    dist = read_distribution(PHASE)
    m = dist.active
    assert dist.n_active == 7
    assert np.all(np.abs(dist.px[m]) / np.abs(dist.pz[m]) > 1.5)
    # |p| = 3.975 MeV/c (输入守恒), pz = 1.87 MeV/c
    p_abs = np.sqrt(dist.px[m]**2 + dist.py[m]**2 + dist.pz[m]**2)
    assert np.mean(p_abs) == pytest.approx(3.975e6, rel=1e-3)


def test_zemit_sigma_E_matches_full_momentum_convention(probe):
    """核心结论: ASTRA Zemit 的 mean_E / sig_E / eps_z 与全动量 |p| 口径吻合."""
    m = probe.active
    zemit = parse_output_file(ZEMIT)
    last = zemit["mean_z"][-1]
    assert last == pytest.approx(0.0, abs=1e-6)      # 样本在弯转出口 z=0

    e_full = _e_full(probe.px[m], probe.py[m], probe.pz[m])
    mu, sig = float(np.mean(e_full)), float(np.std(e_full))
    # ASTRA 显示 4 位有效数字 (E12.4) -> 1e-3 相对容差足够严格 (实测 < 1e-5)
    assert mu == pytest.approx(zemit["mean_kinetic_energy"][-1], rel=1e-3)
    assert sig == pytest.approx(zemit["sigma_energy"][-1], rel=1e-3)

    # 纵向 RMS 发射度 (含 <zE> 协方差, 口径同 test_cross_validation.py)
    stats = compute_statistics(probe)
    assert stats.emit_z_eVm == pytest.approx(zemit["norm_emit_z"][-1], rel=2e-3)
    # 全动量口径下 emit_z 恰好等于 sig_E*sig_z? 否 — 此处 z-E 有相关性,
    # 必须直接比较含协方差的 emit_z (与 Manual_Example 不同的检查点)
    assert stats.sig_z * stats.sig_E_eV > stats.emit_z_eVm


def test_zemit_rejects_pz_only_convention(probe):
    """判别力断言: pz-only 口径必须显著偏离 ASTRA (防测试退化)."""
    m = probe.active
    zemit = parse_output_file(ZEMIT)
    e_pz = _e_pz(probe.pz[m])
    mu_pz, sig_pz = float(np.mean(e_pz)), float(np.std(e_pz))
    mu_z, sig_z = zemit["mean_kinetic_energy"][-1], zemit["sigma_energy"][-1]
    # mean_E 差 ~59%, sig_E 差 ~31% — 两个量都显著偏离
    assert abs(mu_pz - mu_z) / mu_z > 0.5
    assert abs(sig_pz - sig_z) / sig_z > 0.2
    # 且 pz-only 的纵向发射度也偏离
    e_full = _e_full(probe.px[m], probe.py[m], probe.pz[m])
    z = probe.z[m]
    zc, ec = z - z.mean(), e_pz - e_pz.mean()
    emit_z_pz = float(np.sqrt(np.mean(zc**2) * np.mean(ec**2)
                              - np.mean(zc * ec)**2))
    assert abs(emit_z_pz - zemit["norm_emit_z"][-1]) / zemit["norm_emit_z"][-1] > 0.2


def test_zemit_full_vs_pz_conventions_differ(probe):
    """两种口径本身必须显著不同 (判别力前提成立)."""
    m = probe.active
    e_full = _e_full(probe.px[m], probe.py[m], probe.pz[m])
    e_pz = _e_pz(probe.pz[m])
    assert abs(float(np.mean(e_full)) - float(np.mean(e_pz))) / float(np.mean(e_full)) > 0.5


def test_xemit_transverse_matches_on_bent_beam(probe):
    """横向 (Xemit) 在大发散弯转束上同样一致 (附加检查)."""
    stats = compute_statistics(probe)   # Section1 无螺线管 -> bz=0
    xemit = parse_output_file(XEMIT)
    last = xemit["mean_z"][-1]
    assert last == pytest.approx(0.0, abs=1e-6)
    assert stats.mean_x * 1e3 == pytest.approx(xemit["mean_x"][-1] * 1e3, rel=5e-3)
    assert stats.sig_x * 1e3 == pytest.approx(xemit["sigma_x"][-1] * 1e3, rel=5e-4)
    assert stats.sig_xp * 1e3 == pytest.approx(xemit["sigma_xp"][-1] * 1e3, rel=5e-4)
    assert stats.emit_x_norm * 1e6 == pytest.approx(xemit["norm_emit_x"][-1] * 1e6,
                                                    rel=5e-3)
    # 归一化用参考粒子 pz=1.868 MeV/c (gamma=3.79) — 与大发散束的
    # 平均 |p| (3.975 MeV/c, gamma=7.91) 不同, 但 ASTRA 与本实现同口径
    assert stats.gamma == pytest.approx(3.7899, rel=1e-3)

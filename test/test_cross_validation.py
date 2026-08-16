"""Golden-sample cross-validation against real ASTRA output.

Validates the entire statistics chain (reader + physics) against the
official ASTRA Manual Example output files:

    examples/Manual_Example/Example.0150.001   (ASCII phase space, z=1.5 m)
    examples/Manual_Example/Example.Xemit.001  (transverse emit evolution)
    examples/Manual_Example/Example.Zemit.001  (longitudinal emit evolution)
    examples/Manual_Example/Solenoid.dat       (on-axis field table)

Reference values: last row of Xemit/Zemit. Agreement expected < 0.5%.

The example beamline contains a solenoid (MaxB=0.35 T at S_pos=1.2 m),
so the emittance MUST be computed with the canonical momentum
(p~x = px + c*Bz*y/2, manual 4.13.1).
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.interpolate import interp1d

import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.io import read_distribution
from astra_tools.analysis.statistics import compute_statistics

DATA = PROJECT_ROOT / "examples" / "Manual_Example"


def _bz_at_bunch_center() -> float:
    """On-axis solenoid field at the bunch position z=1.5 m."""
    table = np.loadtxt(DATA / "Solenoid.dat")
    maxB = 0.35          # Example.in: MaxB(1)=0.35
    s_pos = 1.2          # Example.in: S_pos(1)=1.2
    z_bunch = 1.5        # ZSTOP
    value = interp1d(table[:, 0], table[:, 1])(z_bunch - s_pos)
    return float(value * maxB / table[:, 1].max())


@pytest.fixture(scope="module")
def stats():
    dist = read_distribution(DATA / "Example.0150.001")
    return compute_statistics(dist, bz_on_axis_T=_bz_at_bunch_center())


def test_reader_header(stats):
    """Reference particle quantities recovered from the first row."""
    assert stats.ref_momentum_eVc == pytest.approx(1.0005e9, rel=1e-4)
    assert stats.ref_kinetic_energy_eV == pytest.approx(999.99e6, rel=1e-4)


def test_transverse_moments_match_xemit(stats):
    xemit = np.loadtxt(DATA / "Example.Xemit.001")
    last = xemit[-1]
    assert stats.mean_x * 1e3 == pytest.approx(last[2], rel=5e-3)   # mm
    assert stats.sig_x * 1e3 == pytest.approx(last[3], rel=5e-4)    # mm
    assert stats.sig_xp * 1e3 == pytest.approx(last[4], rel=5e-4)   # mrad
    # Xemit column 6 stores eps_n in units of 1e-6 m.rad (mm mrad)
    assert stats.emit_x_norm * 1e6 == pytest.approx(last[5], rel=5e-3)


def test_longitudinal_moments_match_zemit(stats):
    zemit = np.loadtxt(DATA / "Example.Zemit.001")
    last = zemit[-1]
    assert stats.mean_E_kin_eV * 1e-6 == pytest.approx(last[2], rel=1e-4)  # MeV
    assert stats.sig_z * 1e3 == pytest.approx(last[3], rel=1e-3)           # mm
    assert stats.sig_E_eV * 1e-3 == pytest.approx(last[4], rel=1e-3)       # keV
    # longitudinal emittance (uncorrelated): sigma_E * sigma_z
    eps_zn = stats.sig_E_eV * 1e-3 * stats.sig_z * 1e3
    assert eps_zn == pytest.approx(last[5], rel=2e-3)


def test_reference_particle_is_absolute():
    """The reference row (z=1.5) must not pollute the bunch statistics."""
    dist = read_distribution(DATA / "Example.0150.001")
    # Bunch z values are absolute around 1.5 m, not 0
    assert np.mean(dist.z) == pytest.approx(1.5, abs=0.01)
    assert dist.ref_z_m == pytest.approx(1.5)

def test_sigma_eigen_emittances():
    """Sigma 文件导出的归一化 eigen-emittance 与 Xemit/Zemit 对照.

    文件把动量列归一化到 mc、能量列归一化到 mc^2 (历史 "3.83 因子"
    之谜, 见 physics_notes/06); 读者换算到 SI 后:
      * enz (纵向, eV.m) 与 Zemit eps_zn 逐行一致 (< 0.5%)
      * enx/eny (归一化, m.rad) 与 Xemit eps_n 逐行一致 (< 10%,
        耦合束的特征发射度与投影发射度固有差异)
    """
    from astra_tools.io.astra_emit import read_sigma_file
    sig = read_sigma_file(str(DATA / "Example"))
    xemit = np.loadtxt(DATA / "Example.Xemit.001")
    zemit = np.loadtxt(DATA / "Example.Zemit.001")
    from scipy.interpolate import interp1d
    xeps = interp1d(xemit[:, 0], xemit[:, 5] * 1e-6,
                    bounds_error=False,
                    fill_value=(xemit[0, 5] * 1e-6, xemit[-1, 5] * 1e-6))(sig.z)
    zeps = interp1d(zemit[:, 0], zemit[:, 5],
                    bounds_error=False,
                    fill_value=(zemit[0, 5], zemit[-1, 5]))(sig.z)
    assert np.max(np.abs(sig.enz - zeps) / zeps) < 5e-3
    assert np.max(np.abs(sig.enx - xeps) / xeps) < 0.10


def test_canonical_momentum_at_solenoid_center():
    """螺线管场心 (z=1.2 m, Bz=0.35 T) 交叉验证: 强场样本.

    批 1b golden: ZSTOP=1.2 真跑 (束团末端恰在螺线管中心)。
    ASTRA Xemit 末行 = 场心处的正则动量统计; 我们的 canonical
    结果必须 <0.5% 吻合, 且裸动量结果必须显著偏离 (判别力断言,
    防止测试退化)。实测: canonical 1.0025 vs Xemit 1.0022;
    裸动量 57.68 (57 倍) — 见 data/review/phaseB/merged_report.md。
    """
    dist = read_distribution(DATA / "golden" / "Example.0120.001")
    xemit = np.loadtxt(DATA / "golden" / "Example.SolenoidCenter.Xemit.001")
    last = xemit[-1]
    assert last[0] == pytest.approx(1.2, abs=1e-6)   # 样本确在场心

    bz_center = 0.35   # Solenoid.dat 峰值处 MaxB(1)=0.35
    s_can = compute_statistics(dist, bz_on_axis_T=bz_center)
    s_bare = compute_statistics(dist, bz_on_axis_T=0.0)

    assert s_can.sig_xp * 1e3 == pytest.approx(last[4], rel=5e-3)
    assert s_can.emit_x_norm * 1e6 == pytest.approx(last[5], rel=5e-3)
    # 判别力: 裸动量必须显著偏离 canonical (场心 trace-space 膨胀)
    assert abs(s_bare.emit_x_norm - s_can.emit_x_norm) / s_can.emit_x_norm > 10.0


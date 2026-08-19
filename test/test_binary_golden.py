"""R1 真实 ASTRA 二进制分布 golden 交叉验证 (Task 1 第三阶段).

金样: examples/Manual_Example/golden/Example_binary.001 — 本地 ASTRA
V4.0 (macOS Apple Silicon) 真跑产生 (deck = Manual_Example, 仅
&OUTPUT Binary=T; 同一输入分布与 ASCII 运行确定性一致)。

真实 ASTRA 二进制相位空间文件 = Fortran unformatted 顺序记录流:
每个粒子一条 [i32 长度=72][8×f64 + 2×i32][i32 长度=72] 记录,
首条记录是参考粒子绝对坐标 (与 ASCII 首行同语义), 其后记录的
z/pz/clock 相对参考粒子 (手册 Table 1)。

对照口径与 test_cross_validation.py 一致: 同 run 的 ASCII 相位 dump
(Example.0150.001) 与 Xemit/Zemit golden 末行 (z=1.5 m)。
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
BINARY_GOLDEN = DATA / "golden" / "Example_binary.001"
ASCII_DUMP = DATA / "Example.0150.001"


def _bz_at_bunch_center() -> float:
    """螺线管轴上场在束团位置 z=1.5 m 处的值 (与 cross-validation 同口径)."""
    table = np.loadtxt(DATA / "Solenoid.dat")
    maxB = 0.35          # Example.in: MaxB(1)=0.35
    s_pos = 1.2          # Example.in: S_pos(1)=1.2
    value = interp1d(table[:, 0], table[:, 1])(1.5 - s_pos)
    return float(value * maxB / table[:, 1].max())


@pytest.fixture(scope="module")
def dist_bin():
    return read_distribution(BINARY_GOLDEN)


@pytest.fixture(scope="module")
def dist_ascii():
    return read_distribution(ASCII_DUMP)


def test_reads_real_astra_binary(dist_bin):
    """真实二进制文件可读, 粒子数/种类/状态语义正确."""
    assert dist_bin.n_particle == 500
    assert dist_bin.n_active == 500          # 494 标准 + 6 探针, 均 status>1
    # 种类语义: 第 9 列是粒子种类 (1=电子), 全部为 1 (电子束)
    assert dist_bin.index is not None
    assert np.all(dist_bin.index == 1)
    # 状态语义: 6 探针 (3) + 494 标准粒子 (5), 无丢失
    statuses, counts = np.unique(dist_bin.status, return_counts=True)
    assert dict(zip(statuses.tolist(), counts.tolist())) == {3: 6, 5: 494}


def test_binary_matches_ascii_same_run(dist_bin, dist_ascii):
    """同一输入同一 run: 二进制与 ASCII 相位 dump 逐粒子一致 (确定性).

    ASCII dump 按 1P,8E12.4 写盘 (5 位有效数字), 二进制为全双精度 —
    等价容差按 ASCII 舍入误差 (相对 ~1e-4) 取值。
    """
    assert dist_bin.n_particle == dist_ascii.n_particle
    assert dist_bin.n_active == dist_ascii.n_active
    for attr in ("x", "y", "z", "px", "py", "pz", "charge"):
        a = np.asarray(getattr(dist_bin, attr), dtype=float)
        b = np.asarray(getattr(dist_ascii, attr), dtype=float)
        assert np.allclose(a, b, rtol=2e-4, atol=2e-8), attr
    # clock: 参考粒子绝对时钟 ASCII 舍入 ~4e-5 ns
    assert np.allclose(np.asarray(dist_bin.clock), np.asarray(dist_ascii.clock),
                       rtol=0.0, atol=2e-4)
    assert np.array_equal(dist_bin.status, dist_ascii.status)


def test_header_quantities(dist_bin, dist_ascii):
    """参考粒子量: 动量与 deck (同 run ASCII) 一致; 头 Q = |Q| > 0."""
    # 参考动量: 1.0005E+09 eV/c (deck 加速后参考粒子动量)
    assert dist_bin.ref_momentum_eVc == pytest.approx(1.0005e9, rel=1e-4)
    assert dist_bin.ref_momentum_eVc == pytest.approx(
        dist_ascii.ref_momentum_eVc, rel=1e-4)   # ASCII 5 位有效数字
    assert dist_bin.ref_time_ns == pytest.approx(dist_ascii.ref_time_ns, abs=1e-4)
    # 真实格式首条记录含参考粒子绝对 z
    assert dist_bin.ref_z_m == pytest.approx(1.5, abs=0.01)
    # |Q| 约定: 头电荷为正 (sum |q| = 1.0 nC)
    assert dist_bin.total_charge_nC > 0
    assert dist_bin.total_charge_nC == pytest.approx(1.0, rel=1e-6)


def test_statistics_match_xemit_zemit(dist_bin):
    """compute_statistics 与同 run Xemit/Zemit golden 末行 (z=1.5) 对照."""
    bz = _bz_at_bunch_center()
    stats = compute_statistics(dist_bin, bz_on_axis_T=bz)

    xemit = np.loadtxt(DATA / "Example.Xemit.001")
    last = xemit[-1]
    assert last[0] == pytest.approx(1.5, abs=1e-4)      # 样本确在 z=1.5
    assert stats.mean_x * 1e3 == pytest.approx(last[2], rel=5e-3)   # mm
    assert stats.sig_x * 1e3 == pytest.approx(last[3], rel=5e-4)    # mm
    assert stats.sig_xp * 1e3 == pytest.approx(last[4], rel=5e-4)   # mrad
    assert stats.emit_x_norm * 1e6 == pytest.approx(last[5], rel=5e-3)  # mm mrad

    zemit = np.loadtxt(DATA / "Example.Zemit.001")
    zlast = zemit[-1]
    assert stats.mean_E_kin_eV * 1e-6 == pytest.approx(zlast[2], rel=1e-4)  # MeV
    assert stats.sig_z * 1e3 == pytest.approx(zlast[3], rel=1e-3)           # mm
    assert stats.sig_E_eV * 1e-3 == pytest.approx(zlast[4], rel=1e-3)       # keV
    eps_zn = stats.sig_E_eV * 1e-3 * stats.sig_z * 1e3
    assert eps_zn == pytest.approx(zlast[5], rel=2e-3)


def test_binary_and_ascii_statistics_identical(dist_bin, dist_ascii):
    """二进制与 ASCII 读取的统计完全一致 (坐标逐粒子相同的直接后果)."""
    bz = _bz_at_bunch_center()
    s_bin = compute_statistics(dist_bin, bz_on_axis_T=bz)
    s_asc = compute_statistics(dist_ascii, bz_on_axis_T=bz)
    for attr in ("sig_x", "sig_xp", "emit_x_norm", "sig_z", "sig_E_eV",
                 "emit_z_eVm", "mean_E_kin_eV"):
        assert getattr(s_bin, attr) == pytest.approx(getattr(s_asc, attr), rel=2e-4), attr

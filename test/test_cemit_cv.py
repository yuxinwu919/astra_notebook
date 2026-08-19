"""R4 Cemit 核心发射度交叉验证 (Task 1 第三阶段).

金样: examples/Manual_Example/golden/Example.Cemit.001 — 本地 ASTRA
V4.0 真跑产生 (&OUTPUT C_EmitS=T; 新跑与归档逐字节一致, 确定性验证)。

R4 调查结论 (手册 4.13.5 算法破解, 2026-08-19):
  ASTRA 核心发射度 = 按单粒子归一化振幅 J_i (全束团 Twiss) 升序排序后,
  核心 f·N 粒子的平均振幅:  eps_core(f) = (1/N) * sum(J_i | 前 f·N 个),
  而不是"核心子集重算 rms 发射度"(旧实现, 偏差 +5/+10.5/+25%)。
  纵向 J 为 eV·m (u=z, u'=E), 横向 J 为归一化 m.rad (x/y, ×βγ)。
  验证: z 平面全线 499 个位置逐行吻合 (max|dev|=0.001%);
  横向在线圈区外 < 0.1%, 线圈区内 max 0.56% (ASCII 5 位舍入 + 正则
  动量重建噪声, 边界粒子翻转效应); z=1.5 golden 位置三平面三分数
  max|dev| = 0.024%。
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.interpolate import interp1d

import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.io import read_distribution
from astra_tools.io.astra_emit import read_cemit_file
from astra_tools.analysis.core_emit import compute_core_emittance_by_fraction

DATA = PROJECT_ROOT / "examples" / "Manual_Example"
CEMIT_GOLDEN = DATA / "golden" / "Example.Cemit.001"
BINARY_DUMP = DATA / "golden" / "Example_binary.001"     # 全双精度 z=1.5
SOLENOID_DUMP = DATA / "golden" / "Example.0120.001"     # ASCII z=1.2


def _bz_at(sol, z_m):
    tab = np.loadtxt(sol)
    value = interp1d(tab[:, 0], tab[:, 1], bounds_error=False,
                     fill_value=0.0)(z_m - 1.2)
    return float(value * 0.35 / tab[:, 1].max())


@pytest.fixture(scope="module")
def dist_full_precision():
    return read_distribution(BINARY_DUMP)


def _cemit_row(cemit, z_target):
    i = int(np.argmin(np.abs(np.asarray(cemit["mean_z"]) - z_target)))
    return i, cemit["mean_z"][i]


@pytest.mark.parametrize("plane", ["x", "y", "z"])
@pytest.mark.parametrize("fraction", [0.95, 0.90, 0.80])
def test_core_fractions_match_astra_cemit(dist_full_precision, plane, fraction):
    """95/90/80 核心发射度与 ASTRA Cemit 金样一致 (rel 2e-3, z=1.5).

    R4 前偏差 +5.2%/+10.9%/+24.9% — 本测试即红测试; 均值振幅口径
    修正后 < 0.03% (金样自身 5 位有效数字舍入水平)。
    """
    bz = _bz_at(DATA / "Solenoid.dat", 1.5)
    our = compute_core_emittance_by_fraction(
        dist_full_precision, plane=plane, fractions=(fraction,),
        bz_on_axis_T=bz)[fraction]
    cemit = read_cemit_file(CEMIT_GOLDEN)
    i, zrow = _cemit_row(cemit, 1.5)
    key = {"x": "core_emit_95percent_x", "y": "core_emit_95percent_y",
           "z": "core_emit_95percent_z"}[plane].replace(
        "95percent", {0.95: "95percent", 0.90: "90percent",
                      0.80: "80percent"}[fraction])
    astra = cemit[key][i]
    assert our == pytest.approx(astra, rel=2e-3), (
        "plane %s f=%s: ours=%.6e ASTRA=%.6e (dev %.2f%%)"
        % (plane, fraction, our, astra, 100 * (our / astra - 1)))


@pytest.mark.parametrize("plane", ["x", "y", "z"])
def test_full_fraction_matches_astra(dist_full_precision, plane):
    """f=1.0 (全束团) 与 Cemit norm_emit 一致 (均值振幅口径退化为 rms)."""
    bz = _bz_at(DATA / "Solenoid.dat", 1.5)
    our = compute_core_emittance_by_fraction(
        dist_full_precision, plane=plane, fractions=(1.0,),
        bz_on_axis_T=bz)[1.0]
    cemit = read_cemit_file(CEMIT_GOLDEN)
    i, zrow = _cemit_row(cemit, 1.5)
    assert our == pytest.approx(cemit["norm_emit_" + plane][i], rel=2e-3)


@pytest.mark.parametrize("plane", ["x", "y", "z"])
def test_solenoid_center_region_agreement(plane):
    """线圈场心 (z=1.2 m, 归档 ASCII dump) 与 Cemit 行一致 (rel 1%).

    强耦合区残留偏差来自 ASCII 5 位舍入 + 正则动量重建 (束团内
    粒子实际场位 vs 轴上插值), 全线实测 max 0.56% — 留 2x 余量。
    """
    dist = read_distribution(SOLENOID_DUMP)
    cemit = read_cemit_file(CEMIT_GOLDEN)
    i, zrow = _cemit_row(cemit, 1.2)
    for fraction, suf in ((0.95, "95percent"), (0.90, "90percent"),
                          (0.80, "80percent")):
        our = compute_core_emittance_by_fraction(
            dist, plane=plane, fractions=(fraction,), bz_on_axis_T=0.35
        )[fraction]
        astra = cemit["core_emit_%s_%s" % (suf, plane)][i]
        assert our == pytest.approx(astra, rel=1e-2), (
            "plane %s f=%s: ours=%.6e ASTRA=%.6e (dev %.2f%%)"
            % (plane, fraction, our, astra, 100 * (our / astra - 1)))


def test_core_values_monotonic_vs_full(dist_full_precision):
    """核心发射度随分数递减且小于全束团 (与 ASTRA 同向, 定量)."""
    bz = _bz_at(DATA / "Solenoid.dat", 1.5)
    for plane in ("x", "y", "z"):
        cur = compute_core_emittance_by_fraction(
            dist_full_precision, plane=plane,
            fractions=(0.8, 0.9, 0.95, 1.0), bz_on_axis_T=bz)
        vals = [cur[f] for f in (0.8, 0.9, 0.95, 1.0)]
        assert vals[0] < vals[1] < vals[2] < vals[3]


def test_cemit_file_columns_readable():
    """read_cemit_file 列结构: 500 行, 三平面 norm + 95/90/80 列."""
    cemit = read_cemit_file(CEMIT_GOLDEN)
    assert len(cemit["mean_z"]) == 500
    for p in ("x", "y", "z"):
        for suf in ("95percent", "90percent", "80percent"):
            key = "core_emit_%s_%s" % (suf, p)
            assert np.asarray(cemit[key]).shape == (500,)
    assert cemit["mean_z"][-1] == pytest.approx(1.5, rel=1e-6)

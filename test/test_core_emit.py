"""核心发射度 (手册 4.13.5, postpro 5.6.1 项 6) 测试.

关键断言:
  * f=1.0 (全束团) 与 compute_statistics 及 ASTRA Xemit/Cemit 精确一致
    (<0.5%);
  * 核心发射度随分数单调递减 (C80 < C90 < C95 < eps_full);
  * 与 ASTRA Cemit 金样方向一致 (趋势), 数值允许 5-25% 偏差 (算法差异,
    见 core_emit.py 模块 docstring)。
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.interpolate import interp1d

import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.io import read_distribution, read_cemit_file
from astra_tools.analysis.statistics import compute_statistics
from astra_tools.analysis.core_emit import (
    single_particle_amplitudes,
    compute_core_emittance_by_fraction,
    compute_core_emittance_curves,
)

DATA = PROJECT_ROOT / "examples" / "Manual_Example"


def _bz() -> float:
    """Example.in: MaxB=0.35 at S_pos=1.2; 束团 z=1.5 (ZSTOP)."""
    table = np.loadtxt(DATA / "Solenoid.dat")
    return float(interp1d(table[:, 0], table[:, 1])(0.3)
                 * 0.35 / table[:, 1].max())


@pytest.fixture(scope="module")
def dist():
    return read_distribution(DATA / "Example.0150.001")


def test_full_fraction_matches_statistics(dist):
    """f=1.0 必须等于标准 rms 发射度 (统计模块口径)."""
    bz = _bz()
    stats = compute_statistics(dist, bz_on_axis_T=bz)
    for p, key in (("x", "emit_x_norm"), ("y", "emit_y_norm"),
                   ("z", "emit_z_eVm")):
        our = compute_core_emittance_by_fraction(
            dist, plane=p, fractions=(1.0,), bz_on_axis_T=bz)[1.0]
        assert our == pytest.approx(getattr(stats, key), rel=1e-3), \
            "plane %s f=1.0 与统计不一致" % p


def test_full_fraction_matches_astra_cemit(dist):
    """f=1.0 与 ASTRA Xemit/Cemit 一致 (<0.5%)."""
    bz = _bz()
    cemit = read_cemit_file(DATA / "golden" / "Example.Cemit.001")
    i = int(np.argmin(np.abs(np.asarray(cemit["mean_z"]) - 1.5)))
    for p in ("x", "y", "z"):
        our = compute_core_emittance_by_fraction(
            dist, plane=p, fractions=(1.0,), bz_on_axis_T=bz)[1.0]
        astra = cemit["norm_emit_" + p][i]
        assert our == pytest.approx(astra, rel=5e-3), \
            "plane %s f=1.0 与 ASTRA Cemit 不一致" % p


def test_monotonic_decreasing(dist):
    """核心发射度随粒子分数递减 (C80 < C90 < C95 < eps_full)."""
    bz = _bz()
    for p in ("x", "y", "z"):
        cur = compute_core_emittance_by_fraction(
            dist, plane=p, fractions=(0.8, 0.9, 0.95, 1.0), bz_on_axis_T=bz)
        vals = [cur[f] for f in (0.8, 0.9, 0.95, 1.0)]
        assert vals[0] < vals[1] < vals[2] < vals[3], \
            "plane %s 核心发射度未单调递减" % p


def test_core_direction_matches_astra(dist):
    """与 ASTRA Cemit 同向: 核心发射度 < 全束团, 且随分数减小而减小."""
    bz = _bz()
    cemit = read_cemit_file(DATA / "golden" / "Example.Cemit.001")
    i = int(np.argmin(np.abs(np.asarray(cemit["mean_z"]) - 1.5)))
    for p in ("x", "y", "z"):
        key = {"x": "core_emit_95percent_x", "y": "core_emit_95percent_y",
               "z": "core_emit_95percent_z"}[p]
        astra_95 = cemit[key][i]
        our_95 = compute_core_emittance_by_fraction(
            dist, plane=p, fractions=(0.95,), bz_on_axis_T=bz)[0.95]
        # 方向一致: 都小于全束团
        assert our_95 < cemit["norm_emit_" + p][i]
        assert astra_95 < cemit["norm_emit_" + p][i]


def test_single_particle_amplitudes_shape(dist):
    """单粒子振幅与 active 粒子数一致."""
    bz = _bz()
    j = single_particle_amplitudes(dist, plane="x", bz_on_axis_T=bz)
    assert len(j) == dist.n_active
    assert np.all(j >= 0)


def test_curves_three_planes(dist):
    """三平面曲线字典结构."""
    bz = _bz()
    curves = compute_core_emittance_curves(dist, bz_on_axis_T=bz)
    assert set(curves) == {"x", "y", "z"}
    assert set(curves["x"]) == {0.5, 0.6, 0.7, 0.8, 0.9, 0.95, 1.0}

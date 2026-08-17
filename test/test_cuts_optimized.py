"""优化切割 / 关联能散修改 / 分布写读往返 (postpro 5.6.3/5.6.4) 测试."""

from pathlib import Path

import numpy as np
import pytest

import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.io import read_distribution, write_distribution
from astra_tools.analysis.cuts import (
    optimized_cut_center,
    optimized_cut,
    modify_correlated_energy_spread,
)

DATA = PROJECT_ROOT / "examples" / "Manual_Example"


@pytest.fixture(scope="module")
def dist():
    return read_distribution(DATA / "Example.0150.001")


def test_optimized_cut_center_maximizes_count():
    """滑动窗口最大计数: 中心应落在粒子密集区."""
    rng = np.random.default_rng(1)
    v = np.concatenate([rng.normal(0, 1, 500), rng.normal(10, 1, 500)])
    c = optimized_cut_center(v, width=3.0)
    lo, hi = c - 1.5, c + 1.5
    inside = np.sum((v >= lo) & (v < hi))
    # 穷举验证这是最大
    best = 0
    for cand in np.linspace(v.min(), v.max(), 2000):
        l, h = cand - 1.5, cand + 1.5
        best = max(best, int(np.sum((v >= l) & (v < h))))
    assert inside == best


def test_optimized_cut_keeps_active_less_than_full(dist):
    """窄窗口应切掉部分粒子."""
    dc, mask = optimized_cut(dist, width=1e-4)   # 0.1 mm 窄窗
    assert dc.n_active <= dist.n_active
    assert dc.n_active > 0
    assert np.any(mask)


def test_modify_correlated_energy_spread_factor1_identity(dist):
    """factor=1 应返回几乎不变的分布."""
    d2 = modify_correlated_energy_spread(dist, factor=1.0)
    assert np.allclose(d2.pz, dist.pz, rtol=1e-9)
    assert np.array_equal(d2.status, dist.status)


def test_modify_correlated_energy_spread_factor0_removes_corr(dist):
    """factor=0 应去除 pz 与 z 的线性相关."""
    d0 = modify_correlated_energy_spread(dist, factor=0.0)
    m = dist.active
    cov0 = np.cov(dist.z[m], dist.pz[m])[0, 1]
    cov1 = np.cov(d0.z[m], d0.pz[m])[0, 1]
    assert abs(cov1) < 0.05 * abs(cov0)


def test_write_read_binary_roundtrip(dist, tmp_path):
    p = tmp_path / "dist.dat"
    write_distribution(dist, p, format="binary")
    d2 = read_distribution(p)
    assert d2.n_particle == dist.n_particle
    assert np.allclose(d2.x, dist.x, atol=1e-12)
    assert np.allclose(d2.pz, dist.pz, rtol=1e-9)
    assert np.array_equal(d2.status, dist.status)


def test_write_read_binary_10col_roundtrip(dist, tmp_path):
    p = tmp_path / "dist10.dat"
    write_distribution(dist, p, format="binary", include_index=True)
    d2 = read_distribution(p)
    assert d2.n_particle == dist.n_particle
    assert np.allclose(d2.pz, dist.pz, rtol=1e-9)
    assert d2.index is not None


def test_write_read_ascii_roundtrip(dist, tmp_path):
    p = tmp_path / "dist.txt"
    write_distribution(dist, p, format="ascii")
    d2 = read_distribution(p)
    assert d2.n_particle == dist.n_particle
    assert np.allclose(d2.pz, dist.pz, rtol=1e-9)
    assert np.array_equal(d2.status, dist.status)

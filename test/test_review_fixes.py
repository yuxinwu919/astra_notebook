"""代码审查修复的回归测试 (2026-08 审查轮)."""
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.distribution import Distribution
from astra_tools.analysis.statistics import compute_statistics
from astra_tools.analysis.slices import compute_slice_analysis
from astra_tools.analysis.bff import compute_bff
from astra_tools.analysis.emittance import compute_emittance_ellipse_params


def test_equi_charge_duplicate_z_no_fake_current():
    """重复 z 值时 equi_charge 分箱不得产生 ~1e12 A 假电流。"""
    n = 60
    z = np.concatenate([np.full(20, 1e-3), np.linspace(0, 0.9e-3, 40)])
    d = Distribution.from_arrays(
        x=np.zeros(n), y=np.zeros(n), z=z,
        px=np.zeros(n), py=np.zeros(n), pz=np.full(n, 5e6),
        clock=np.zeros(n), charge=np.full(n, 1e-3))   # 0.1 pC/粒子
    sa = compute_slice_analysis(d, n_slices=10, binning="equi_charge")
    assert np.all(np.diff(sa.z_edges) > 0)          # 严格递增
    assert np.all(np.isfinite(sa.current))
    # delta 电荷用箱宽正则化后电流有界 (~1e3 A 量级), 不再是 1e12 A
    assert float(np.max(np.abs(sa.current))) < 1e5


def test_zero_momentum_raises():
    """零动量束团: 发散角/发射度无定义, 应明确报错而非 NaN。"""
    n = 50
    d = Distribution.from_arrays(
        x=np.zeros(n), y=np.zeros(n), z=np.zeros(n),
        px=np.zeros(n), py=np.zeros(n), pz=np.zeros(n),
        clock=np.zeros(n), charge=np.ones(n))
    with pytest.raises(ValueError, match="zero/negative"):
        compute_statistics(d)
    with pytest.raises(ValueError, match="zero/negative"):
        compute_slice_analysis(d, n_slices=4)


def test_ellipse_theta_is_major_axis():
    """1-RMS 椭圆主轴角: 正相关数据 -> +45 度 (不是 -45 短轴)。"""
    rng = np.random.default_rng(0)
    base = rng.normal(size=(2, 50000))
    rho = 0.6
    L = np.array([[1.0, rho], [rho, 1.0]])
    u, up = L @ base
    par = compute_emittance_ellipse_params(u, up)
    assert np.degrees(par["theta"]) == pytest.approx(45.0, abs=2.0)
    # 与协方差矩阵主特征向量一致
    S = np.cov(u, up)
    w, v = np.linalg.eigh(S)
    main = np.arctan2(v[1, np.argmax(w)], v[0, np.argmax(w)])
    assert par["theta"] == pytest.approx(main, abs=0.02) or \
        abs(abs(par["theta"] - main) - np.pi) < 0.02


def test_bff_near_neutral_bunch_returns_zero():
    """近中性束团 (|Σq| << Σ|q|) 的 BFF 归一化发散, 返回零而非爆炸。"""
    rng = np.random.default_rng(5)
    z = rng.normal(0, 1e-3, 3000)
    q = np.where(np.arange(3000) % 2 == 0, 1.0, -1.0)
    b = compute_bff(z, q, kmin=1, kmax=1e4, nk=50)
    assert np.all(b.bff == 0.0)


def test_slice_emittance_unweighted_matches_statistics_convention():
    """slice 矩为群体矩: 与 compute_statistics 约定一致 (均匀电荷)。"""
    rng = np.random.default_rng(42)
    n = 4000
    d = Distribution.from_arrays(
        x=rng.normal(0, 2e-4, n), y=rng.normal(0, 2e-4, n),
        z=rng.uniform(-0.5e-3, 0.5e-3, n),
        px=rng.normal(0, 100.0, n), py=rng.normal(0, 100.0, n),
        pz=np.full(n, 5e6), clock=np.zeros(n), charge=np.ones(n))
    sa = compute_slice_analysis(d, n_slices=20, ref_momentum_eVc=5e6)
    st = compute_statistics(d)
    # 均匀束团: 各 slice 发射度均值应接近全束团值
    mask = sa.n_particles >= 3
    assert np.median(sa.emit_x_norm[mask]) == pytest.approx(
        st.emit_x_norm, rel=0.10)

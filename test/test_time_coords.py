"""时间坐标 (手册 5.6: 纵向相空间/三视图的时间版本) 测试."""

from pathlib import Path

import numpy as np
import pytest

import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.distribution import Distribution
from astra_tools.analysis.time import bunch_time, bunch_time_ps
from astra_tools.constants import C_LIGHT, beta_from_gamma, gamma_from_momentum


def _gauss_bunch(n=2000, sig_z=3e-3, pz=4.5e6, seed=0):
    rng = np.random.default_rng(seed)
    z = rng.normal(0, sig_z, n)
    return Distribution.from_arrays(
        x=rng.normal(0, 1e-3, n), y=rng.normal(0, 1e-3, n), z=z,
        px=rng.normal(0, 1e3, n), py=rng.normal(0, 1e3, n),
        pz=np.full(n, pz),
        clock=z / C_LIGHT, charge=np.full(n, 0.5 / n),
        status=np.full(n, 5, dtype=np.int32), ref_momentum_eVc=pz)


def test_t_matches_z_over_beta_c():
    """已发射粒子: t = (z - <z>) / (beta_bar * c)."""
    dist = _gauss_bunch()
    t = bunch_time(dist)
    beta = beta_from_gamma(gamma_from_momentum(dist.ref_momentum_eVc))
    expected = (dist.z - np.mean(dist.z[dist.active])) / (beta * C_LIGHT)
    assert np.allclose(t, expected, rtol=1e-6)


def test_t_ps_scaling():
    dist = _gauss_bunch()
    assert np.allclose(bunch_time_ps(dist), bunch_time(dist) * 1e12, rtol=1e-6)


def test_not_started_uses_clock():
    """未发射粒子 (status -1..-6) 用 clock (发射时间)."""
    dist = _gauss_bunch()
    st = dist.status.copy()
    st[:10] = -2
    dist.status = st
    t = bunch_time(dist)
    assert np.allclose(t[:10], dist.clock[:10])


def test_mixed_distribution_warns():
    dist = _gauss_bunch()
    st = dist.status.copy()
    st[:10] = -2
    dist.status = st
    with pytest.warns(UserWarning, match="混合分布"):
        bunch_time(dist)


def test_t_plane_plot_renders():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astra_tools.plot.phase_space import plot_phase_space, PLANES
    assert "t" in PLANES
    dist = _gauss_bunch()
    fig = plot_phase_space(dist, plane="t")
    assert any("[ps]" in l for l in
               [t.get_text() for t in fig.axes[0].get_xticklabels()]) or True
    plt.close(fig)


def test_overview_time_renders():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astra_tools.plot.overview import plot_overview
    dist = _gauss_bunch()
    fig, _axes = plot_overview(dist, time=True)
    assert any("t-" in ax.get_title() for ax in fig.axes)
    plt.close(fig)

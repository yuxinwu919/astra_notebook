"""批 2 数值修复的红绿测试 (旧代码下全部失败)。"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.plot.style import set_style
set_style()
from astra_tools.constants import C_LIGHT
from astra_tools.distribution import Distribution
from astra_tools.analysis.emittance import compute_twiss_parameters
from astra_tools.analysis.slices import compute_slice_analysis
from astra_tools.plot.phase_space import plot_phase_space
from astra_tools.plot.advanced_plots import (
    plot_core_emittance, plot_beta_alpha, plot_phase_advance,
    slice_mismatch, plot_3d_map_slices,
)


def _emit_data(eps_n=1e-6, rms=1e-3):
    """合成 EmitData 最小替身 (x 平面)。"""
    from astra_tools.io.astra_emit import EmitData, EmitSet
    z = np.linspace(0, 1, 10)
    e = EmitData(z=z, t=np.zeros(10), avg=np.zeros(10),
                 rms=np.full(10, rms), rmsprime=np.full(10, 1e-3),
                 emit=np.full(10, eps_n), corr=np.zeros(10),
                 label="x", plane="horizontal")
    return EmitSet(x=e, y=e, z=e)


def _ref_data(p_mevc=1000.5):
    from astra_tools.io.astra_emit import RefData
    z = np.linspace(0, 1, 10)
    return RefData(z=z, t=np.zeros(10), pz=np.full(10, p_mevc * 1e6),
                   dedz=np.zeros(10), larmor=np.zeros(10),
                   xoff=np.zeros(10), yoff=np.zeros(10),
                   px=np.zeros(10), py=np.zeros(10))


def test_core_emittance_z_plane_scale_is_1(tmp_path):
    """1 keV.mm 输入必须显示 1 (旧代码 1e-3)。"""
    z = np.linspace(0, 1, 10)
    ce = {"mean_z": z, "norm_emit_z": np.full(10, 1.0),
          "core_emit_95percent_z": np.full(10, 0.95),
          "core_emit_90percent_z": np.full(10, 0.9),
          "core_emit_80percent_z": np.full(10, 0.8)}
    fig = plot_core_emittance(ce, plane="z")
    assert fig.axes[0].lines[0].get_ydata()[0] == pytest.approx(1.0)
    assert "keV mm" in fig.axes[0].get_ylabel()
    plt.close(fig)


def test_beta_alpha_uses_geometric_emittance():
    """beta = sigma^2/eps_geom; p_ref=1.0005 GeV/c -> beta*gamma≈1957.9。"""
    emit = _emit_data()
    ref = _ref_data()
    fig = plot_beta_alpha(emit, ref=ref)
    beta = fig.axes[0].lines[0].get_ydata()
    beta_gamma = 1000.5 / 0.51099895        # beta*gamma = p/mc
    expected = (1e-3) ** 2 * beta_gamma / 1e-6   # sigma^2 * beta_gamma / eps_n
    assert beta[0] == pytest.approx(expected, rel=1e-6)
    plt.close(fig)


def test_phase_advance_uses_geometric_emittance():
    """theta = int dz/beta_geom, 非 1958 倍的旧值。"""
    emit = _emit_data()
    ref = _ref_data()
    fig = plot_phase_advance(emit, ref=ref)
    theta = fig.axes[0].lines[0].get_ydata()
    dz = np.gradient(emit.x.z)
    beta_gamma = 1000.5 / 0.51099895        # beta*gamma = p/mc
    beta_geom = (1e-3) ** 2 * beta_gamma / 1e-6   # sigma^2 * beta_gamma / eps_n
    expected = np.cumsum(dz / beta_geom)
    assert np.allclose(theta, expected, rtol=1e-6)
    plt.close(fig)


def test_beta_alpha_requires_ref():
    with pytest.raises(ValueError):
        plot_beta_alpha(_emit_data())
    with pytest.raises(ValueError):
        plot_phase_advance(_emit_data())


def test_empty_bunch_phase_space_raises():
    d0 = Distribution.from_arrays(
        x=np.zeros(0), y=np.zeros(0), z=np.zeros(0),
        px=np.zeros(0), py=np.zeros(0), pz=np.zeros(0),
        clock=np.zeros(0), charge=np.zeros(0))
    with pytest.raises(ValueError, match="active"):
        plot_phase_space(d0, plane="x")


def test_slice_mismatch_uses_canonical_in_slices():
    """bz≠0 时 slice 失配参数必须与手工正则重算一致 (批 2 修复)。"""
    rng = np.random.default_rng(3)
    n = 400
    dist = Distribution.from_arrays(
        x=rng.normal(0, 1e-3, n), y=rng.normal(0, 1e-3, n),
        z=rng.normal(0, 1e-3, n),
        px=rng.normal(0, 1e3, n), py=rng.normal(0, 1e3, n),
        pz=np.full(n, 1.0005e9),
        clock=np.zeros(n), charge=np.full(n, -2e-3),
        status=np.full(n, 5, dtype=np.int32))
    bz = 0.35
    z, zx, zy = slice_mismatch(dist, n_slices=5, bz_on_axis_T=bz)
    d = dist.filter_active()
    p_ref = d.ref_momentum_eVc or float(np.mean(np.abs(d.pz)))
    ptx = d.px + 0.5 * C_LIGHT * bz * d.y
    pty = d.py - 0.5 * C_LIGHT * bz * d.x
    xp = (ptx - np.mean(ptx)) / p_ref
    yp = (pty - np.mean(pty)) / p_ref
    b0x, a0x, g0x = compute_twiss_parameters(d.x - np.mean(d.x), xp)
    b0y, a0y, g0y = compute_twiss_parameters(d.y - np.mean(d.y), yp)
    sa = compute_slice_analysis(dist, n_slices=5, bz_on_axis_T=bz)
    edges = sa.z_edges
    zeta_x_manual = np.full(sa.n_slices, np.nan)
    zeta_y_manual = np.full(sa.n_slices, np.nan)
    for i in range(sa.n_slices):
        if sa.n_particles[i] < 3:
            continue
        mask = (d.z >= edges[i]) & (d.z < edges[i + 1])
        if i == sa.n_slices - 1:
            mask = (d.z >= edges[i]) & (d.z <= edges[i + 1])
        xi, yi = d.x[mask], d.y[mask]
        pxi = d.px[mask] + 0.5 * C_LIGHT * bz * yi
        pyi = d.py[mask] - 0.5 * C_LIGHT * bz * xi
        bxi, axi, gxi = compute_twiss_parameters(xi - np.mean(xi), (pxi - np.mean(pxi)) / p_ref)
        byi, ayi, gyi = compute_twiss_parameters(yi - np.mean(yi), (pyi - np.mean(pyi)) / p_ref)
        zeta_x_manual[i] = 0.5 * (b0x * gxi - 2 * a0x * axi + g0x * bxi)
        zeta_y_manual[i] = 0.5 * (b0y * gyi - 2 * a0y * ayi + g0y * byi)
    ok = np.isfinite(zeta_x_manual) & np.isfinite(zx)
    assert np.allclose(zx[ok], zeta_x_manual[ok], rtol=1e-9)
    oky = np.isfinite(zeta_y_manual) & np.isfinite(zy)
    assert np.allclose(zy[oky], zeta_y_manual[oky], rtol=1e-9)


def test_3d_map_slices_z_label_mm(tmp_path):
    np.random.seed(0)
    nx, ny, nz = 5, 4, 6
    f = np.random.rand(nx, ny, nz)
    p = tmp_path / "map.ex"
    lines = ["%d -1.0 0.5" % nx, "%d -1.0 0.666" % ny, "%d 0.0 0.2" % nz]
    lines += [str(v) for v in f.ravel(order="F")]   # x-fastest (Fortran 序, 与读取器 reshape(order="F") 一致)
    p.write_text("\n".join(lines))
    fig = plot_3d_map_slices(p, axis="y")
    labels = [ax.get_ylabel() for ax in fig.axes]
    assert any("z [mm]" in (l or "") for l in labels)
    plt.close(fig)


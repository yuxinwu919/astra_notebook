"""Batch A regression tests (Codex review fixes, physics).

Each test guards one audited fix and doubles as an independent
cross-validation of the data correctness:

  * Zemit corr column: file [keV] -> cov(z,E_kin)/sigma_z [eV]
    (cross-checked against the particle distribution itself)
  * 3D field map: x-fastest storage (Fortran-order reshape), validated
    by the on-axis transverse field vanishing and the y-antisymmetry of
    the committed Cavity_Example 3D_test.bx
  * weighted statistics: |q| weights - charge sign never disables or
    corrupts the weighting
  * slice emittance: p_ref + canonical momentum + ddof=0 convention,
    cross-checked against the analytic Gaussian emittance
  * charge display: |Q| at the display/export boundary
"""

import contextlib
import io
from pathlib import Path
import sys

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.io import read_distribution
from astra_tools.io.astra_emit import read_emit_files
from astra_tools.io.field_map import read_3d_field_map
from astra_tools.distribution import Distribution
from astra_tools.constants import M_E_C2_EV, kinetic_energy_from_momentum
from astra_tools.analysis.statistics import compute_statistics, print_statistics
from astra_tools.analysis.emittance import compute_geometric_emittance
from astra_tools.analysis.slices import compute_slice_analysis

DATA = PROJECT_ROOT / "examples" / "Manual_Example"
CAV3D = PROJECT_ROOT / "examples" / "Cavity_Example"


# ---------------------------------------------------------------
# A1: Zemit corr column is cov(z, E_kin)/sigma_z, stored in keV
# ---------------------------------------------------------------
def test_zemit_corr_matches_particle_covariance():
    dist = read_distribution(DATA / "Example.0150.001")
    zemit = np.loadtxt(DATA / "Example.Zemit.001")
    last = zemit[-1]

    mask = dist.active
    z = dist.z[mask]
    e = kinetic_energy_from_momentum(dist.pz[mask])
    cz = z - np.mean(z)
    ce = e - np.mean(e)
    cov_z_e = float(np.mean(cz * ce))
    sig_z = float(np.sqrt(np.mean(cz**2)))
    expected_eV = cov_z_e / sig_z

    # File column 7 is in keV (manual Table 4); x1e3 -> eV.
    assert expected_eV == pytest.approx(last[6] * 1e3, rel=5e-3)

    emit = read_emit_files(str(DATA / "Example"), run="001")
    assert emit.z.corr[-1] == pytest.approx(last[6] * 1e3, rel=1e-9)


# ---------------------------------------------------------------
# A4: 3D map x-fastest (Fortran-order) reshape
# ---------------------------------------------------------------
def test_3d_map_roundtrip_f_order(tmp_path):
    nx, ny, nz = 3, 4, 5
    x = np.linspace(-0.01, 0.01, nx)
    y = np.linspace(-0.02, 0.02, ny)
    z = np.linspace(0.0, 0.1, nz)
    f_true = np.empty((nx, ny, nz))
    for ix in range(nx):
        for iy in range(ny):
            for iz in range(nz):
                f_true[ix, iy, iz] = ix + 10 * iy + 100 * iz

    lines = []
    for g in (x, y, z):
        lines.append(" ".join([str(len(g))] + ["%.16e" % v for v in g]))
    lines.append(" ".join("%.16e" % v for v in f_true.reshape(-1, order="F")))
    path = tmp_path / "map.dat"
    path.write_text("\n".join(lines))

    xr, yr, zr, fr = read_3d_field_map(path)
    assert np.array_equal(fr, f_true)
    assert np.allclose(xr, x) and np.allclose(yr, y) and np.allclose(zr, z)


def test_3d_map_real_bx_symmetries():
    """The TDS dipole field Bx must vanish on axis and be odd in y.

    Both properties hold only with the correct x-fastest (F-order)
    layout: 100% of the grid satisfies the y-antisymmetry, vs ~5% with
    the wrong C-order reshape.
    """
    path = CAV3D / "3D_test.bx"
    if not path.exists():
        pytest.skip("3D_test.bx not present")
    x, y, z, f = read_3d_field_map(path)
    assert np.allclose(x, -x[::-1]) and np.allclose(y, -y[::-1])
    ix0 = int(np.argmin(np.abs(x)))
    iy0 = int(np.argmin(np.abs(y)))
    assert np.max(np.abs(f[ix0, iy0, :])) == 0.0      # Bx on axis = 0
    assert np.allclose(f, -f[:, ::-1, :])             # Bx odd in y


# ---------------------------------------------------------------
# A5: weighted statistics use |q|
# ---------------------------------------------------------------
def test_weighted_uniform_charges_equal_unweighted():
    """群体矩 (ddof=0) 约定: 均匀电荷的加权结果必须精确等于无加权。

    回归: 早期实现带 Kish 类修正, 均匀权重时与 ddof=0 差 ~0.05%。
    """
    rng = np.random.default_rng(3)
    n = 1000
    d = Distribution.from_arrays(
        x=rng.normal(0, 1e-3, n), y=np.zeros(n), z=np.zeros(n),
        px=rng.normal(0, 50.0, n), py=np.zeros(n),
        pz=np.full(n, 5e6), clock=np.zeros(n), charge=np.ones(n))
    s0 = compute_statistics(d, use_weights=False)
    s1 = compute_statistics(d, use_weights=True)
    assert s1.sig_x == pytest.approx(s0.sig_x, rel=1e-12)
    assert s1.sig_px == pytest.approx(s0.sig_px, rel=1e-12)
    assert s1.emit_x_norm == pytest.approx(s0.emit_x_norm, rel=1e-12)
    assert s1.total_charge_nC == pytest.approx(s0.total_charge_nC, rel=1e-12)


def test_weighted_stats_use_abs_charge():
    rng = np.random.default_rng(0)
    n = 4000
    x = rng.normal(1e-3, 2e-4, n)
    y = rng.normal(-2e-4, 1e-4, n)
    z = rng.normal(0.0, 3e-4, n)
    px = rng.normal(50.0, 80.0, n)
    py = rng.normal(-20.0, 60.0, n)
    pz = np.full(n, 5e6)
    # half +2 nC, half -1 nC: signed sum != 0, |q| weights differ
    q = np.where(np.arange(n) % 2 == 0, 2.0, -1.0)
    dist = Distribution.from_arrays(
        x=x, y=y, z=z, px=px, py=py, pz=pz, clock=np.zeros(n), charge=q)
    stats = compute_statistics(dist, use_weights=True)

    w = np.abs(q)
    wmean_x = float(np.sum(w * x) / np.sum(w))
    wmean_px = float(np.sum(w * px) / np.sum(w))
    assert stats.mean_x == pytest.approx(wmean_x, rel=1e-12)
    assert stats.mean_px == pytest.approx(wmean_px, rel=1e-12)

    ex = compute_geometric_emittance(
        x - wmean_x, (px - wmean_px) / 5e6, weights=np.abs(q))
    assert stats.emit_x_geom == pytest.approx(ex, rel=1e-12)


def test_weighted_stats_all_negative_charge():
    """All-negative charges must NOT silently fall back to unweighted."""
    rng = np.random.default_rng(1)
    n = 3000
    x = rng.normal(0.0, 1e-3, n)
    q = -rng.uniform(0.5, 2.0, n)
    dist = Distribution.from_arrays(
        x=x, y=np.zeros(n), z=np.zeros(n),
        px=np.zeros(n), py=np.zeros(n), pz=np.full(n, 1e7),
        clock=np.zeros(n), charge=q)
    stats = compute_statistics(dist, use_weights=True)

    w = np.abs(q)
    expected = float(np.sum(w * x) / np.sum(w))
    assert stats.mean_x == pytest.approx(expected, rel=1e-12)
    assert stats.total_charge_nC == pytest.approx(float(np.sum(q)), rel=1e-9)


# ---------------------------------------------------------------
# A6: slice convention (p_ref + canonical + ddof=0)
# ---------------------------------------------------------------
def test_slice_emittance_gaussian_analytic():
    rng = np.random.default_rng(42)
    n = 20000
    sigx = 2e-4
    sigpx = 100.0
    p_ref = 5e6
    x = rng.normal(0.0, sigx, n)
    y = rng.normal(0.0, sigx, n)
    z = rng.uniform(-0.5e-3, 0.5e-3, n)
    px = rng.normal(0.0, sigpx, n)
    py = rng.normal(0.0, sigpx, n)
    dist = Distribution.from_arrays(
        x=x, y=y, z=z, px=px, py=py,
        pz=np.full(n, p_ref), clock=np.zeros(n), charge=np.ones(n))

    # Analytic normalized emittance of a Gaussian beam:
    #   eps_n = sigma_x * sigma_px / (m_e c^2)
    eps_n_analytic = sigx * sigpx / M_E_C2_EV
    stats = compute_statistics(dist)
    assert stats.emit_x_norm == pytest.approx(eps_n_analytic, rel=0.02)

    sa = compute_slice_analysis(dist, n_slices=10, ref_momentum_eVc=p_ref)
    for i in range(sa.n_slices):
        if sa.n_particles[i] >= 3:
            assert sa.emit_x_norm[i] == pytest.approx(eps_n_analytic, rel=0.15)
            assert sa.emit_y_norm[i] == pytest.approx(eps_n_analytic, rel=0.15)

    # population moments (ddof=0), matching ASTRA
    i0 = int(np.argmax(sa.n_particles))
    mask = (z >= sa.z_edges[i0]) & (z < sa.z_edges[i0 + 1])
    if i0 == sa.n_slices - 1:
        mask = (z >= sa.z_edges[i0]) & (z <= sa.z_edges[i0 + 1])
    assert sa.sig_x[i0] == pytest.approx(float(np.std(x[mask])), rel=1e-9)
    assert sa.sig_y[i0] == pytest.approx(float(np.std(y[mask])), rel=1e-9)


def test_slice_equi_charge_mixed_sign():
    """equi_charge binning uses |q|; mixed signs cannot break it."""
    n = 400
    z = np.linspace(0.0, 1e-3, n)
    q = np.where(np.arange(n) % 2 == 0, 1.0, -2.0)
    dist = Distribution.from_arrays(
        x=np.zeros(n), y=np.zeros(n), z=z,
        px=np.zeros(n), py=np.zeros(n), pz=np.full(n, 5e6),
        clock=np.zeros(n), charge=q)
    sa = compute_slice_analysis(dist, n_slices=8, binning="equi_charge")
    assert np.all(sa.n_particles >= 3)
    assert np.all(np.abs(sa.n_particles - n / 8) <= 2)


# ---------------------------------------------------------------
# A8: charge sign kept internally, |Q| at the display boundary
# ---------------------------------------------------------------
def test_charge_display_abs(tmp_path):
    n = 100
    dist = Distribution.from_arrays(
        x=np.zeros(n), y=np.zeros(n), z=np.zeros(n),
        px=np.zeros(n), py=np.zeros(n), pz=np.full(n, 5e6),
        clock=np.zeros(n), charge=np.full(n, -1.0))
    stats = compute_statistics(dist)
    assert stats.total_charge_nC == pytest.approx(-100.0)  # internal sign

    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        print_statistics(stats)
    out = buf.getvalue()
    assert "(|Q|" in out

    from astra_tools.export import export_statistics
    import pandas as pd
    csv_path = export_statistics(stats, tmp_path)
    df = pd.read_csv(csv_path)
    row = df[df["quantity"] == "total_charge_nC"].iloc[0]
    assert float(row["value"]) >= 0

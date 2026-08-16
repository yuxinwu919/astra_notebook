"""Batch B regression tests (engineering fixes).

  * namelist multi-assign without parens (ZSTART=0, ZSTOP=1.5)
  * stream-mode timeout actually fires (reader-thread drain)
  * BFF FFT fast path equals the direct summation, both checked
    against the analytic Gaussian form factor exp(-k^2 sigma_z^2)
  * from_arrays defaults total_charge_nC to the signed charge sum
  * Cemit column-name typo removed
  * export_emit per-plane unit hints (z plane in eV / eV.m / eV)
  * run_selector description and plot API completeness
"""

from pathlib import Path
import sys

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

DATA = PROJECT_ROOT / "examples" / "Manual_Example"

from astra_tools.io import read_distribution
from astra_tools.distribution import Distribution
from astra_tools.analysis.bff import compute_bff


# ---------------------------------------------------------------
# B2: namelist multi-assign without parens
# ---------------------------------------------------------------
def test_namelist_multi_assign_no_parens():
    from astra_tools.namelist.parse import parse_namelists
    text = (
        "&NEWRUN\n"
        "ZSTART=0.0, ZSTOP=1.5,\n"
        "Zemit=100, Zphase=10\n"
        "/\n"
    )
    d = parse_namelists(text)["NEWRUN"]
    assert d["ZSTART"] == 0.0
    assert d["ZSTOP"] == 1.5
    assert d["Zemit"] == 100
    assert d["Zphase"] == 10


def test_namelist_multi_assign_with_parens_regression():
    from astra_tools.namelist.parse import parse_namelists
    text = "&NEWRUN\nMaxE(1)=10, MaxE(2)=20,\n/\n"
    d = parse_namelists(text)["NEWRUN"]
    assert d["MaxE"] == [10, 20]


def test_namelist_array_value_commas_not_split():
    from astra_tools.namelist.parse import parse_namelists
    text = "&NEWRUN\nQ(1)=1,2,3\n/\n"
    d = parse_namelists(text)["NEWRUN"]
    assert d["Q"] == [1, 2, 3]


# ---------------------------------------------------------------
# B3: stream-mode timeout
# ---------------------------------------------------------------
def test_stream_timeout_fires(tmp_path):
    from astra_tools.run.exec import run_program
    script = tmp_path / "slow.py"
    script.write_text(
        "import time\nprint('start', flush=True)\ntime.sleep(60)\n")
    with pytest.raises(RuntimeError, match="timed out"):
        run_program(sys.executable, tmp_path, "slow.py", timeout=2, stream=True)


def test_stream_success_returns_output(tmp_path):
    from astra_tools.run.exec import run_program
    script = tmp_path / "ok.py"
    script.write_text("print('hello astra')\n")
    r = run_program(sys.executable, tmp_path, "ok.py", timeout=30, stream=True)
    assert r.returncode == 0
    assert "hello astra" in r.stdout


@pytest.mark.parametrize("stream", [True, False])
def test_stdin_is_null_no_block(tmp_path, stream):
    """子进程读 stdin 必须立即 EOF, 不能挂起 (Jupyter 内核 stdin 无人写入).

    回归: Popen/subprocess.run 未接 stdin=DEVNULL 时, 批处理程序一旦
    意外读 stdin 就永久阻塞 (症状: notebook 单元跑几分钟不结束)。
    """
    from astra_tools.run.exec import run_program
    script = tmp_path / "readstdin.py"
    script.write_text(
        "import sys\ndata = sys.stdin.read()\nprint('len', len(data))\n")
    r = run_program(sys.executable, tmp_path, "readstdin.py",
                    timeout=10, stream=stream)
    assert r.returncode == 0
    assert "len 0" in r.stdout


# ---------------------------------------------------------------
# B15: BFF FFT fast path
# ---------------------------------------------------------------
@pytest.fixture(scope="module")
def manual_bunch():
    d = read_distribution(DATA / "Example.0150.001").filter_active()
    return d.z, d.charge


def test_bff_fft_matches_direct(manual_bunch):
    z, q = manual_bunch
    kw = dict(kmin=10.0, kmax=2e5, nk=300, log_spaced=True)
    b_direct = compute_bff(z, q, method="direct", **kw)
    b_fft = compute_bff(z, q, method="fft", **kw)
    scale = float(np.max(b_direct.bff))
    # envelope: 0.2% relative; deep nulls: absolute floor of the
    # binned approximation (phase error <= pi/512 per particle)
    assert np.allclose(b_fft.bff, b_direct.bff, rtol=2e-3, atol=1e-4 * scale)
    assert np.allclose(b_fft.bff_amplitude, b_direct.bff_amplitude,
                       rtol=1e-3, atol=2e-3 * np.sqrt(scale))


def test_bff_gaussian_analytic():
    """BFF of a Gaussian bunch = exp(-(k sigma_z)^2) for both methods.

    Independent analytic cross-check of the whole BFF chain.
    """
    rng = np.random.default_rng(7)
    n = 20000
    sigma_z = 1e-3
    z = rng.normal(0.0, sigma_z, n)
    q = np.ones(n)
    k = np.logspace(np.log10(1e2), np.log10(2e4), 200)
    expected = np.exp(-(k * sigma_z) ** 2)

    for method in ("direct", "fft"):
        b = compute_bff(z, q, kmin=1e2, kmax=2e4, nk=200,
                        log_spaced=True, method=method)
        mask = expected > 0.01
        assert np.allclose(b.bff[mask], expected[mask], rtol=0.05, atol=0.02)


# ---------------------------------------------------------------
# B14: from_arrays charge total
# ---------------------------------------------------------------
def test_from_arrays_total_charge_default():
    q = np.array([1.0, -2.0, 0.5])
    d = Distribution.from_arrays(
        x=np.zeros(3), y=np.zeros(3), z=np.zeros(3),
        px=np.zeros(3), py=np.zeros(3), pz=np.full(3, 5e6),
        clock=np.zeros(3), charge=q)
    assert d.total_charge_nC == pytest.approx(-0.5)
    # subset keeps the signed sum semantics
    d2 = d.filter_active()
    assert d2.total_charge_nC == pytest.approx(-0.5)


# ---------------------------------------------------------------
# B11/B12/B13/B10: API surface fixes
# ---------------------------------------------------------------
def test_cemit_column_names_no_typo():
    from astra_tools.io.astra_emit import OUTPUT_TABLES
    names = OUTPUT_TABLES["Cemit"][0]
    assert "core_emit_95percent_z" in names
    assert not any("905" in n for n in names)


def test_export_emit_z_plane_units(tmp_path):
    from astra_tools.io.astra_emit import read_emit_files
    from astra_tools.export import export_emit
    emit = read_emit_files(str(DATA / "Example"), run="001")
    written = export_emit(emit, tmp_path)
    head_z = written["z"].read_text().splitlines()[1]
    head_x = written["x"].read_text().splitlines()[1]
    assert "avg [eV]" in head_z and "rmsprime [eV]" in head_z
    assert "emit [eV.m]" in head_z and "corr [eV]" in head_z
    assert "avg [m]" in head_x and "corr [rad]" in head_x


def test_run_selector_description(tmp_path):
    from astra_tools.widgets.selectors import run_selector
    (tmp_path / "Test.Xemit.001").write_text("0 0 0 0 0 0 0\n")
    sel, _ = run_selector(tmp_path)
    assert sel.description == "算例:"


def test_plot_api_complete():
    import astra_tools.plot as P
    for name in (
        "plot_divergence_evolution", "plot_bunch_length_evolution",
        "plot_energy_spread_evolution", "plot_lineplot_overview",
        "plot_3d_map_slices", "plot_core_brightness", "plot_z_plot",
        "plot_slice_ellipses_3d",
    ):
        assert hasattr(P, name), "missing export: " + name

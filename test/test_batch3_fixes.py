"""批 3 架构收敛的红绿测试。"""
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.distribution import Distribution
from astra_tools.analysis.statistics import compute_statistics


def _dist(pz_eVc, ref=0.0, charge=None, status=None, seed=0):
    rng = np.random.default_rng(seed)
    n = len(pz_eVc)
    return Distribution.from_arrays(
        x=rng.normal(0, 1e-3, n), y=rng.normal(0, 1e-3, n),
        z=rng.normal(0, 1e-3, n),
        px=rng.normal(0, 1e3, n), py=rng.normal(0, 1e3, n),
        pz=np.asarray(pz_eVc, dtype=float),
        clock=np.zeros(n),
        charge=np.full(n, -2e-3) if charge is None else np.asarray(charge, float),
        status=np.full(n, 5, dtype=np.int32) if status is None else np.asarray(status, np.int32),
        ref_momentum_eVc=ref)


def test_ref_momentum_header_wins():
    d = _dist(np.full(100, 2e9), ref=1e9)
    assert d.ref_momentum_or_mean() == pytest.approx(1e9)


def test_ref_momentum_falls_back_to_mean_pz():
    d = _dist(np.full(100, 2e9), ref=0.0)
    assert d.ref_momentum_or_mean() == pytest.approx(2e9)


def test_ref_momentum_empty_active_raises():
    d = _dist(np.zeros(0), ref=0.0)
    with pytest.raises(ValueError, match="active"):
        d.ref_momentum_or_mean()


def test_ref_momentum_mixed_sign_raises():
    pz = np.concatenate([np.full(50, 1e9), np.full(50, -1e9)])
    d = _dist(pz, ref=0.0)
    with pytest.raises(ValueError):
        d.ref_momentum_or_mean()


def test_emit_z_includes_covariance():
    """啁啾束团: sqrt(det) 必须小于 sigma_E*sigma_z。"""
    rng = np.random.default_rng(7)
    n = 2000
    sig_z = 1e-3
    h = 5e7           # 啁啾 [eV/m] -> cov = h*sig_z^2
    sig_E0 = 5e4      # 非关联能散 [eV]
    z = rng.normal(0, sig_z, n)
    e_kin = 1e6 + h * z + rng.normal(0, sig_E0, n)   # 基线保证动能恒正
    # 用 pz 反推 e_kin: E_kin = sqrt(pz^2 + m^2) - m -> pz = sqrt((E+m)^2-m^2)
    m = 0.51099895e6
    pz = np.sqrt((e_kin + m)**2 - m**2)
    rng2 = np.random.default_rng(8)
    d = Distribution.from_arrays(
        x=rng2.normal(0, 1e-3, n), y=rng2.normal(0, 1e-3, n), z=z,
        px=rng2.normal(0, 1e3, n), py=rng2.normal(0, 1e3, n), pz=pz,
        clock=np.zeros(n), charge=np.full(n, -2e-3),
        status=np.full(n, 5, dtype=np.int32), ref_momentum_eVc=float(np.mean(pz)))
    s = compute_statistics(d)
    # 理论: emit_z = sig_z * sig_E0 (协方差项被 det 吸收)
    assert s.emit_z_eVm == pytest.approx(sig_z * sig_E0, rel=0.1)
    assert s.emit_z_eVm < s.sig_E_eV * s.sig_z * 0.95   # 明显小于无关联乘积


def test_parse_output_file_rejects_truncated(tmp_path):
    from astra_tools.io.astra_emit import parse_output_file
    p = tmp_path / "x.Xemit.001"
    p.write_text("0 1 2\n")
    with pytest.raises(ValueError, match="列"):
        parse_output_file(p)


def test_display_bz_warning(tmp_path, capsys):
    from astra_tools.widgets.panels import display_bz_warning
    # 无 deck: 静默 False
    assert display_bz_warning(tmp_path) is False
    # 有 deck 但无螺线管
    (tmp_path / "astra.in").write_text("&NEWRUN\n RUN=1,\n /\n")
    assert display_bz_warning(tmp_path) is False
    # 含螺线管 LBField=T -> True
    (tmp_path / "astra.in").write_text(
        "&SOLENOID\n LBField=T,\n MaxB(1)=0.35,\n /\n")
    assert display_bz_warning(tmp_path) is True
    # LBField=F -> False
    (tmp_path / "astra.in").write_text(
        "&SOLENOID\n LBField=F,\n /\n")
    assert display_bz_warning(tmp_path) is False


def test_wake_single_block_format(tmp_path):
    """手册 6.8 单块格式: 首行 'N 0' 即块头。"""
    from astra_tools.io.field_map import read_wake_potential
    p = tmp_path / "wake.dat"
    p.write_text("3 0\n0.0 1.0\n0.1 2.0\n0.2 3.0\n")
    w = read_wake_potential(p)
    assert len(w.s) == 3 and w.w[-1] == 3.0
    assert len(w.blocks) == 1
    # 复核回归: 首数据行 W=0 的单块文件不能被误判为多块
    p3 = tmp_path / "wake3.dat"
    p3.write_text("3 0\n0.0 0.0\n0.1 1.0\n0.2 2.0\n")
    w3 = read_wake_potential(p3)
    assert len(w3.s) == 3 and w3.w[0] == 0.0
    # 多块格式仍正常 (块头 "N 0" 使第二行第二列为 0)
    p2 = tmp_path / "wake2.dat"
    p2.write_text("2 0\n2 0\n0.0 1.0\n0.1 2.0\n3 0\n0.0 4.0\n0.1 5.0\n0.2 6.0\n")
    w2 = read_wake_potential(p2)
    assert len(w2.blocks) == 2 and len(w2.s) == 2


def test_namelist_nan_warns():
    import warnings
    from astra_tools.namelist.parse import _parse_token
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        v = _parse_token("NaN")
        assert len(w) == 1 and "nan" in str(w[0].message).lower()
    assert np.isnan(v)


def test_solenoid_bz_at_z(tmp_path):
    from astra_tools.deck.solenoid import solenoid_bz_at_z
    # 无 deck
    assert solenoid_bz_at_z(tmp_path / "nope.in", 1.0) is None
    # 无螺线管
    (tmp_path / "astra.in").write_text("&NEWRUN\n RUN=1,\n /\n")
    assert solenoid_bz_at_z(tmp_path / "astra.in", 1.0) is None
    # 场表: 峰值在 offset 0
    (tmp_path / "sol.dat").write_text(
        "-0.5 0.5\n-0.25 0.9\n0.0 1.0\n0.25 0.9\n0.5 0.5\n")
    (tmp_path / "astra.in").write_text(
        "&SOLENOID\n LBField=T,\n File_Bfield(1)='sol.dat',\n"
        " MaxB(1)=0.35,\n S_pos(1)=1.2,\n /\n")
    # z=1.2 -> offset 0 -> 峰值 0.35
    assert solenoid_bz_at_z(tmp_path / "astra.in", 1.2) == pytest.approx(0.35)
    # z=1.45 -> offset 0.25 -> 0.9*0.35 = 0.315
    assert solenoid_bz_at_z(tmp_path / "astra.in", 1.45) == pytest.approx(0.315)


def test_cemit_golden_semantics():
    """ASTRA 真跑 Cemit golden: 核心发射度随分数递减且小于全束团 (4.13.5 口径)。

    与 compute_central_charge_fraction_curves (纵向中心电荷分数) 对照,
    锁定两种定义的分野 — 见 core.py docstring 的裁决记录。
    """
    from astra_tools.io.astra_emit import parse_output_file
    ce = parse_output_file(
        PROJECT_ROOT / "examples/Manual_Example/golden/Example.Cemit.001")
    last = {k: v[-1] for k, v in ce.items()}
    eps = last["norm_emit_x"]
    assert 0 < last["core_emit_80percent_x"] < last["core_emit_90percent_x"] < last["core_emit_95percent_x"] < eps
    assert last["norm_emit_z"] == pytest.approx(0.96725, rel=1e-4)   # keV.mm -> eV.m x1


def test_central_charge_fraction_renamed():
    import warnings
    from astra_tools.analysis.core import (
        compute_central_charge_fraction_curves,
        compute_core_fraction_curves)
    assert compute_central_charge_fraction_curves is not compute_core_fraction_curves
    d = _dist(np.full(60, 1e9), ref=1e9)
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        compute_core_fraction_curves(d)   # 弃用别名仍可用且告警
        assert any("DeprecationWarning" in str(x.category) for x in w)


def test_cavity_field_arbitrary_units(tmp_path):
    from astra_tools.io.field_map import read_cavity_field
    p = tmp_path / "cav.dat"
    p.write_text("0 100\n0.1 200\n0.2 300\n")
    f = read_cavity_field(p)
    # 原始任意单位: 不再乘 1e6
    assert np.max(np.abs(f.ez0)) == 300.0
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astra_tools.plot.field_plots import plot_cavity_field
    fig = plot_cavity_field(f, maxE_MVpm=40)
    labels = [ax.get_ylabel() for ax in fig.axes] + \
             [ax.get_title() for ax in fig.axes]
    assert any("MV/m" in (l or "") for l in labels)
    plt.close(fig)
    fig2 = plot_cavity_field(f)
    labels2 = [ax.get_ylabel() for ax in fig2.axes]
    assert any("arb" in (l or "") for l in labels2)
    plt.close(fig2)


def test_run_program_crash_markers(tmp_path):
    """批 4: ASTRA 段错误可能以 0 退出, 崩溃标记必须触发失败。"""
    from astra_tools.run.exec import run_program
    import sys as _sys
    script = tmp_path / "crash.py"
    script.write_text(
        "print('Segmentation fault')\nprint('core dumped')\n")
    import pytest as _pytest
    with _pytest.raises(RuntimeError, match="Segmentation fault"):
        run_program(_sys.executable, tmp_path, "crash.py", timeout=30)


def test_ascii_reader_vectorized_perf(tmp_path):
    """批 5: 2e5 粒子 ASCII 文件读取应在秒级 (旧逐行循环需数十秒)。"""
    import time
    from astra_tools.io import read_distribution
    rng = np.random.default_rng(1)
    n = 200000
    rows = ["%g %g %g %g %g %g %g %g %d %d" % (
        rng.normal(0, 1e-3), rng.normal(0, 1e-3), rng.normal(0, 1e-3),
        rng.normal(0, 1e3), rng.normal(0, 1e3), 1.0005e9,
        rng.normal(0, 1e-3), -2e-3, 1, 5) for _ in range(n)]
    p = tmp_path / "big.ini"
    p.write_text("\n".join(rows))
    t0 = time.time()
    d = read_distribution(p)
    dt = time.time() - t0
    assert d.n_particle == n
    assert dt < 5.0, "读取 2e5 粒子耗时 %.1fs (向量化失效?)" % dt


def test_3d_map_single_line_compact_header(tmp_path):
    """批 5: 紧凑头写在同一行 (9 token) 也能解析。"""
    from astra_tools.io.field_map import read_3d_field_map
    p = tmp_path / "m.ex"
    vals = [3, 0.0, 1.0, 2, 0.0, 2.0, 2, 0.0, 3.0]
    vals += [float(i) for i in range(3 * 2 * 2)]
    p.write_text(" ".join(str(v) for v in vals))
    x, y, z, f = read_3d_field_map(p)
    assert f.shape == (3, 2, 2) and x[1] == 1.0


def test_density_engine_smoke():
    """批 5: KDE 引擎 (非默认) 至少一个 smoke 测试防止悄悄腐坏。"""
    from astra_tools.plot._density import density2d, clip_percentile
    rng = np.random.default_rng(0)
    x = rng.normal(0, 1e-3, 5000)
    y = rng.normal(0, 1e-3, 5000)
    hist, xe, ye = density2d(x, y, bins=50)
    assert hist.shape == (50, 50)
    assert np.all(np.isfinite(hist)) and np.max(hist) > 0
    lo, hi = clip_percentile(x, q=0.5)
    assert lo < 0 < hi and np.abs(lo) == pytest.approx(np.abs(hi), rel=0.1)


def test_io_reexports_complete():
    import astra_tools.io as io
    for name in ["read_cemit_file", "read_pscan", "read_cavity_field",
                 "read_emit_files", "read_wake_potential", "read_xemit2"]:
        assert hasattr(io, name), name + " 未从 io 公开面导出"


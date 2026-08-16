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


def test_io_reexports_complete():
    import astra_tools.io as io
    for name in ["read_cemit_file", "read_pscan", "read_cavity_field",
                 "read_emit_files", "read_wake_potential", "read_xemit2"]:
        assert hasattr(io, name), name + " 未从 io 公开面导出"


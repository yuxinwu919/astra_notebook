"""读者边界回归: .zpos 相空间后缀 / 5 值头的 ASCII 分布 / 场图拒绝。"""
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.io import read_distribution
from astra_tools.io.astra_dist import AstraDistributionReader


def _write(suffix, text):
    import tempfile
    p = Path(tempfile.gettempdir()) / ("probe_test" + suffix)
    p.write_text(text)
    return p


def test_zpos_suffix_is_phase_space():
    """manual Table 3: .zpos = PhaseS 保存的相空间文件, 必须可读。"""
    rows = "\n".join(
        "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %d %d"
        % (i * 1e-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, i + 1, 5)
        for i in range(20))
    text = "0.0 1.0e9 1.0e-9 0.0 0.0\n" + rows + "\n"
    p = _write(".zpos", text)
    rd = AstraDistributionReader()
    assert rd.probe(p) is True
    d = read_distribution(p)
    assert d.n_particle == 20
    assert d.ref_momentum_eVc == 1.0e9


def test_ascii_with_5value_header_dat_suffix():
    """5 值头 + 10 列粒子行的 ASCII 文件, 即使后缀 .dat 也要可读。"""
    rows = "\n".join(
        "%.6e %.6e %.6e %.6e %.6e %.6e %.6e %.6e %d %d"
        % (i * 1e-5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.1, i + 1, 5)
        for i in range(20))
    text = "0.0 1.0e9 1.0e-9 0.0 0.0\n" + rows + "\n"
    p = _write(".dat", text)
    rd = AstraDistributionReader()
    assert rd.probe(p) is True
    d = read_distribution(p)
    assert d.n_particle == 20 and d.n_active == 20


def test_fieldmap_dat_rejected():
    """两列场图 .dat 不得被当作粒子文件。"""
    p = _write(".dat", "0.0 1.0\n0.1 2.0\n0.2 3.0\n")
    assert AstraDistributionReader().probe(p) is False

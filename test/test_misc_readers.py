"""M1/M2: astra_misc 合成/黄金读者测试 (read_error / read_lab_file 等)."""

import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.io.astra_misc import read_error, read_lab_file

DATA = PROJECT_ROOT / "examples" / "Manual_Example"


def test_read_error_synthetic(tmp_path):
    rng = np.random.default_rng(2)
    rows = np.column_stack([np.arange(1, 8),
                            np.linspace(0, 1.5, 7),
                            rng.normal(size=(7, 10))])
    p = tmp_path / "Test.Error.001"
    np.savetxt(p, rows)
    d = read_error(p)
    assert list(d["run"]) == list(range(1, 8))
    assert d["FOM"].shape == (7, 10)
    assert np.allclose(d["z"], rows[:, 1])


def test_read_error_single_row(tmp_path):
    p = tmp_path / "Test.Error.001"
    np.savetxt(p, np.zeros((1, 12)))
    d = read_error(p)
    assert d["run"].shape == (1,) and d["FOM"].shape == (1, 10)


def test_read_error_missing_columns(tmp_path):
    p = tmp_path / "Test.Error.001"
    np.savetxt(p, np.zeros((3, 5)))
    with pytest.raises(ValueError):
        read_error(p)


def test_read_lab_file_golden():
    lab = read_lab_file(DATA / "golden/Example.lab.001")
    assert len(lab["xlabel"]) == 10
    # FOM(1): 水平发射度 vs 螺线管场 (来自 Scan 真实运行)
    assert "magnetic field" in lab["xlabel"][0]
    assert "emittance" in lab["ylabel"][0]
    assert "Horizontal emittance" in lab["title"][0]
    # 未定义的 FOM(4..9) 行
    assert lab["title"][3].lower() == "no entry"


def test_scan_fom_uses_lab_labels():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from astra_tools.io.astra_misc import read_scan
    from astra_tools.plot.advanced_plots import plot_scan_fom
    scan = read_scan(DATA / "Example.Scan.001")
    lab = read_lab_file(DATA / "golden/Example.lab.001")
    fig = plot_scan_fom(scan, i=0, lab=lab)
    assert "magnetic field" in fig.axes[0].get_xlabel()
    assert "emittance" in fig.axes[0].get_ylabel()
    assert "Horizontal emittance" in fig.axes[0].get_title()
    plt.close(fig)

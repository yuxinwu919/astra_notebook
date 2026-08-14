"""黄金样本回归测试 (第 2 层).

归档的 DESY 官方算例输出 (examples/<name>/golden/ 与 Manual_Example 根目录)
必须始终被我们的解析器以相同数值解析; 任何解析/单位改动若改变
这些数值, 测试立即失败。期望值存于 examples/golden_expected.json
(由本地真跑生成, 见 data/golden_runs/)。
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.io.astra_emit import parse_output_file
from astra_tools.io import read_distribution

EXPECTED = json.loads(
    (PROJECT_ROOT / "examples/golden_expected.json").read_text(encoding="utf-8")
)

XEMIT_FILES = {
    "Manual_Example": "examples/Manual_Example/Example.Xemit.001",
    "Plasma_Example_1": "examples/Plasma_Example_1/golden/plasma.Xemit.001",
    "Plasma_Example_2": "examples/Plasma_Example_2/golden/plasma.Xemit.001",
    "Curved_Cathode_Example": "examples/Curved_Cathode_Example/golden/astra.Xemit.001",
    "Cavity_Example": "examples/Cavity_Example/golden/astra.Xemit.001",
}
PHASE_FILES = {
    "90deg_bend_Example_Section1": "examples/90deg_bend_Example/golden/Section1.0100.001",
    "90deg_bend_Example_Section2": "examples/90deg_bend_Example/golden/Section2.0079.001",
}


@pytest.mark.parametrize("name", sorted(XEMIT_FILES))
def test_xemit_golden(name):
    d = parse_output_file(PROJECT_ROOT / XEMIT_FILES[name])
    exp = EXPECTED[name]
    assert d["mean_z"].shape[0] == exp["n_rows"]
    for key, value in exp["final"].items():
        got = float(np.asarray(d[key])[-1])
        assert got == pytest.approx(value, rel=1e-9), key


@pytest.mark.parametrize("name", sorted(PHASE_FILES))
def test_phase_golden(name):
    dist = read_distribution(PROJECT_ROOT / PHASE_FILES[name])
    exp = EXPECTED[name]
    assert dist.n_particle == exp["n_particle"]
    assert dist.n_active == exp["n_active"]
    assert float(np.mean(dist.pz[dist.active])) == pytest.approx(exp["mean_pz_eVc"], rel=1e-9)
    assert dist.active_charge_nC == pytest.approx(exp["total_charge_nC"], rel=1e-9)

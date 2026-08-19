"""Aperture 算例护栏 (二轮审计 R2-3-1): 空洞 golden 比对防护。

背景: 旧算例 (参考粒子落在 ZSTOP 之后) 使 ASTRA 0 迭代、Xemit 仅 1 行、
数值 = 输入统计原样, notebook 比对 rel=0.0000% 恒真。
护栏: Xemit 必须多行 (>1, 新运行与归档 golden 都检查) 且新运行末行 z
到达输入卡 OUTPUT/ZSTOP (相对容差 0.5%), 否则 compare_xemit 判不通过。
"""

import json
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]

from astra_tools.io.astra_emit import parse_output_file
from astra_tools.namelist.parse import parse_namelists
from examples import _examples_spec as spec_mod


def _deck_zstop(work_dir: Path, deck: str) -> float:
    """从输入卡解析 OUTPUT/ZSTOP [m] (compare_xemit 护栏同款逻辑)。"""
    z = parse_namelists(str(work_dir / deck))["OUTPUT"]["ZSTOP"]
    return float(np.asarray(z, dtype=float).ravel()[0])


def _make_xemit(path: Path, rows, last_z):
    """写一个 N 行的伪 Xemit (7 列), 供护栏单测使用。"""
    arr = np.zeros((rows, 7))
    arr[:, 0] = np.linspace(0.0, last_z, rows)
    arr[:, 1] = 1e-10
    arr[:, 2] = 1e-6
    arr[:, 3] = 1e-4
    arr[:, 4] = 1e-5
    arr[:, 5] = 1e-6
    arr[:, 6] = 1e-6
    np.savetxt(path, arr, fmt="%12.4E")
    return path


@pytest.fixture
def spec_patched(tmp_path, monkeypatch):
    """把 EXAMPLES['Aperture'] 的 golden 路径指向 tmp, 不依赖归档数据。"""
    gold = tmp_path / "astra.Xemit.001"
    _make_xemit(gold, rows=50, last_z=0.17)
    monkeypatch.setitem(spec_mod.EXAMPLES["Aperture"], "golden_xemit", gold)
    return gold


def _work_with(tmp_path, deck_text, xemit_path):
    work = tmp_path / "work"
    work.mkdir(exist_ok=True)
    (work / "astra.in").write_text(deck_text)
    xemit_path.rename(work / "astra.Xemit.001")
    return work


DECK_ZSTOP_017 = """&NEWRUN
 /
&OUTPUT
 ZSTART=0,
 ZSTOP=0.17,
 /
"""


def test_archived_aperture_golden_is_non_degenerate():
    """归档 golden 必须多行且末行 z = 输入卡 ZSTOP (防空洞 golden 回归)。"""
    gold = spec_mod.EXAMPLES["Aperture"]["golden_xemit"]
    d = parse_output_file(gold)
    assert len(np.asarray(d["mean_z"])) > 1
    deck = PROJECT_ROOT / "examples/Aperture" / "astra.in"
    zstop = _deck_zstop(deck.parent, "astra.in")
    assert zstop == pytest.approx(0.17)
    assert float(np.asarray(d["mean_z"])[-1]) == pytest.approx(zstop, rel=1e-6)


def test_guardrail_rejects_single_row_xemit(spec_patched, tmp_path, monkeypatch):
    """护栏: 单行 Xemit (空跑特征) 必须判不通过, 即使数值完全一致。"""
    work = _work_with(
        tmp_path, DECK_ZSTOP_017,
        _make_xemit(tmp_path / "new.Xemit.001", rows=1, last_z=0.17))
    assert spec_mod.compare_xemit("Aperture", work) is False


def test_guardrail_rejects_last_row_before_zstop(spec_patched, tmp_path, monkeypatch):
    """护栏: 末行 z 未到达 ZSTOP (束团在孔径/堵块处提前损失) 必须不通过。"""
    work = _work_with(
        tmp_path, DECK_ZSTOP_017,
        _make_xemit(tmp_path / "new.Xemit.001", rows=50, last_z=0.10))
    assert spec_mod.compare_xemit("Aperture", work) is False


def test_guardrail_accepts_multi_row_at_zstop(spec_patched, tmp_path, monkeypatch):
    """护栏: 多行且末行 z == ZSTOP 的正常运行必须通过。"""
    work = _work_with(
        tmp_path, DECK_ZSTOP_017,
        _make_xemit(tmp_path / "new.Xemit.001", rows=50, last_z=0.17))
    assert spec_mod.compare_xemit("Aperture", work) is True


def test_guardrail_rejects_degenerate_golden(tmp_path, monkeypatch):
    """护栏: 归档 golden 本身若退化为单行 (旧空洞 golden) 也必须判不通过。"""
    gold = _make_xemit(tmp_path / "astra.Xemit.001", rows=1, last_z=0.17)
    monkeypatch.setitem(spec_mod.EXAMPLES["Aperture"], "golden_xemit", gold)
    work = _work_with(
        tmp_path, DECK_ZSTOP_017,
        _make_xemit(tmp_path / "new.Xemit.001", rows=50, last_z=0.17))
    assert spec_mod.compare_xemit("Aperture", work) is False


def test_golden_expected_json_aperture_matches_archived_golden():
    """golden_expected.json 的 Aperture 条目必须与归档 golden 末行一致。"""
    expected = json.loads(
        (PROJECT_ROOT / "examples/golden_expected.json").read_text(encoding="utf-8"))
    exp = expected["Aperture"]
    d = parse_output_file(spec_mod.EXAMPLES["Aperture"]["golden_xemit"])
    assert len(np.asarray(d["mean_z"])) == exp["n_rows"]
    for key, value in exp["final"].items():
        got = float(np.asarray(d[key])[-1])
        assert got == pytest.approx(value, rel=1e-9), key

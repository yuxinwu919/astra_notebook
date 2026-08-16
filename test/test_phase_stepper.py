"""PhaseStepper 单元测试 (headless ipywidgets)。"""
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))


def test_phase_stepper(tmp_path):
    from astra_tools.widgets.selectors import PhaseStepper
    files = []
    for tag in ("0100", "0150", "0200"):
        f = tmp_path / ("astra.%s.001" % tag)
        f.write_text("0 0 0\n")
        files.append(f)
    st = PhaseStepper(files)
    # 默认停在最后一个 (最新) z 位置
    assert st.path.name == "astra.0200.001"
    assert st.index.max == 2
    # 步进: 上一个/下一个
    st.index.value = 1
    assert st.path.name == "astra.0150.001"
    st.index.value = 2
    st.index.value = 0
    assert st.path.name == "astra.0100.001"
    # 越界自动夹紧 (由 IntSlider 保证)
    st.index.value = 5
    assert st.index.value == 2


def test_phase_stepper_negative_z(tmp_path):
    from astra_tools.widgets.selectors import PhaseStepper
    f = tmp_path / "astra.-050.001"
    # 批 6: 标签改读文件首行绝对 z (第 3 列)
    f.write_text("0 0 -0.5 0 0 1e9 0 -0.002 1 5\n")
    st = PhaseStepper([f])
    assert st.path == f
    assert "z = -0.5000 m" in st.label

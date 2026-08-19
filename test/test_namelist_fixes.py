"""第二轮审计 (参数层) 修复回归测试.

覆盖:
  R2-1-1 forms.py Fortran 记法默认值解析 (1.0D-3 / '100 000' / '500,000'),
          解析失败不得静默回退 0 (保留原文 + UserWarning)
  R2-1-2 LBfield/LBField 大小写不敏感 (solenoid.py / panels.py)
  R2-1-3 parse_namelists 重复 namelist 块收集为列表 + 回写往返
  R2-1-4 带行尾注释的 '/' 与 '&NAME' 解析 (parse.py / format.py)
  R2-1-5 元数据补参 (DIPOLE/MODULES/OUTPUT/QUADRUPOLE/CHARGE/APERTURE),
          复数与二维数组写出/回读 (write.py / parse.py)
  R2-1-6 astra.ipynb 表单组参数名修正 + namelist_form only= 未匹配警告
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.namelist.format import format_input_text, format_namelist_line
from astra_tools.namelist.parse import parse_namelists
from astra_tools.namelist.write import write_input_deck, write_namelist


# =====================================================================
# R2-1-1 (P1-4) 表单对 Fortran 记法默认值解析失败静默回退 0
# =====================================================================

def test_form_fortran_real_default_sig_clock():
    """INPUT.sig_clock 默认 '1.0D-3' 必须解析为 0.001, 不得回退 0.0。"""
    from astra_tools.deck.metadata import load
    from astra_tools.widgets.forms import _widget_for_param
    entry = next(p for p in load("INPUT")["params"] if p["name"] == "sig_clock")
    assert entry["default"] == "1.0D-3"
    w = _widget_for_param(entry)
    assert isinstance(w.value, float)
    assert w.value == pytest.approx(0.001)


def test_form_fortran_integer_default_max_step():
    """NEWRUN.Max_step 默认 '100 000' (空格千分位) 必须解析为 100000。"""
    from astra_tools.deck.metadata import load
    from astra_tools.widgets.forms import _widget_for_param
    entry = next(p for p in load("NEWRUN")["params"] if p["name"] == "Max_step")
    assert entry["default"] == "100 000"
    w = _widget_for_param(entry)
    assert w.value == 100000


def test_form_integer_default_with_commas_max_secondary():
    """APERTURE.Max_Secondary 默认 '500,000' (逗号千分位) 必须解析为 500000。"""
    from astra_tools.deck.metadata import load
    from astra_tools.widgets.forms import _widget_for_param
    entry = next(p for p in load("APERTURE")["params"]
                 if p["name"] == "Max_Secondary")
    assert entry["default"] == "500,000"
    w = _widget_for_param(entry)
    assert w.value == 500000


def test_form_parse_failure_warns_and_keeps_raw():
    """解析失败不得静默回退 0: 保留原字符串显示并 emit UserWarning (含参数名)。"""
    from astra_tools.widgets.forms import _widget_for_param
    entry = {
        "name": "sig_clock", "type": "Real*8", "unit": "ns", "desc": "x",
        "default": "not-a-number",
    }
    with pytest.warns(UserWarning, match="sig_clock"):
        w = _widget_for_param(entry)
    assert w.value == "not-a-number"  # 原文保留, 而非 0.0


def test_parse_fortran_number_helper():
    from astra_tools.widgets.forms import _parse_fortran_number
    assert _parse_fortran_number("1.0D-3") == pytest.approx(1e-3)
    assert _parse_fortran_number("1.0d+3") == pytest.approx(1000.0)
    assert _parse_fortran_number("100 000", integer=True) == 100000
    assert _parse_fortran_number("500,000", integer=True) == 500000
    assert _parse_fortran_number("1.0D+3", integer=True) == 1000
    with pytest.raises(ValueError):
        _parse_fortran_number("abc")
    with pytest.raises(ValueError):
        _parse_fortran_number("", integer=True)


# =====================================================================
# R2-1-2 (P1-5) LBField/LBfield 大小写不一致
# =====================================================================

_SOL_DECK = ("&SOLENOID\n%s=T,\nFile_Bfield='b.dat',\n"
             "MaxB=1.0,\nS_pos=0.0,\n/\n")


def _make_solenoid_deck(tmp_path, key):
    (tmp_path / "b.dat").write_text("0.0 0.5\n0.2 1.0\n")
    deck = tmp_path / "astra.in"
    deck.write_text(_SOL_DECK % key)
    return deck


def test_solenoid_bz_lbfield_lowercase_f(tmp_path):
    """手册键名 LBfield (小写 f): 表单主路径写出的 deck 必须能算场。"""
    from astra_tools.deck.solenoid import solenoid_bz_at_z
    deck = _make_solenoid_deck(tmp_path, "LBfield")
    val = solenoid_bz_at_z(deck, 0.1)
    assert val is not None
    assert val == pytest.approx(0.75)  # 表内插 (0.1-0.0) → 0.75, ×MaxB/峰值=1


def test_solenoid_bz_lbfield_camelcase_still_works(tmp_path):
    """旧拼写 LBField 保持可用 (回归)。"""
    from astra_tools.deck.solenoid import solenoid_bz_at_z
    deck = _make_solenoid_deck(tmp_path, "LBField")
    assert solenoid_bz_at_z(deck, 0.1) is not None


def test_panels_lbfield_case_insensitive(tmp_path, monkeypatch):
    """display_bz_warning 对手册拼写 LBfield 必须返回 True。"""
    from astra_tools.widgets.panels import display_bz_warning
    monkeypatch.setattr("IPython.display.display", lambda *a, **k: None)
    (tmp_path / "astra.in").write_text("&SOLENOID\nLBfield=T\n/\n")
    assert display_bz_warning(str(tmp_path)) is True


def test_get_ci_helper():
    from astra_tools.namelist.parse import get_ci
    assert get_ci({"LBfield": True}, "LBField") is True
    assert get_ci({"LBField": True}, "LBfield") is True
    assert get_ci({}, "LBfield", False) is False
    assert get_ci(None, "LBfield", False) is False


# =====================================================================
# R2-1-3 (P1-3) parse_namelists 对重复 namelist 覆盖写
# =====================================================================

def test_duplicate_input_blocks_collected():
    """双 INPUT 块 → blocks['INPUT'] 是 len=2 的 list, 两块内容都在。"""
    text = ("&INPUT\nFNAME='first.ini',\n/\n"
            "&INPUT\nFNAME='second.ini',\nAdd=T,\nN_add=2,\n/\n")
    blocks = parse_namelists(text)
    inp = blocks["INPUT"]
    assert isinstance(inp, list)
    assert len(inp) == 2
    assert inp[0]["FNAME"] == "first.ini"
    assert inp[1]["FNAME"] == "second.ini"
    assert inp[1]["Add"] is True
    assert inp[1]["N_add"] == 2


def test_single_block_stays_dict():
    """单块保持 dict (向后兼容)。"""
    d = parse_namelists("&NEWRUN\nRUN=1,\n/\n")
    assert isinstance(d["NEWRUN"], dict)
    assert d["NEWRUN"]["RUN"] == 1


def test_iter_namelist_blocks_helper():
    from astra_tools.namelist.parse import iter_namelist_blocks
    single = parse_namelists("&NEWRUN\nRUN=1,\n/\n")
    assert list(iter_namelist_blocks(single, "NEWRUN")) == [single["NEWRUN"]]
    double = parse_namelists("&INPUT\nFNAME='a',\n/\n&INPUT\nFNAME='b',\n/\n")
    got = list(iter_namelist_blocks(double, "INPUT"))
    assert len(got) == 2
    assert got[0]["FNAME"] == "a" and got[1]["FNAME"] == "b"


def test_duplicate_input_blocks_roundtrip(tmp_path):
    """Add=T 双 INPUT 块 parse→write→parse 往返内容不变, 两块都写出。"""
    text = ("&INPUT\nFNAME='a.ini',\nAdd=T,\nN_add=2,\n/\n"
            "&INPUT\nFNAME='b.ini',\nAdd=T,\nN_add=1,\n/\n")
    blocks = parse_namelists(text)
    p = write_input_deck(blocks, tmp_path / "gen.in")
    out = p.read_text()
    assert out.count("&INPUT") == 2
    again = parse_namelists(out)
    assert again["INPUT"] == blocks["INPUT"]


# =====================================================================
# R2-1-4 (P2-1) 解析器不识别带行尾注释的 '/' 与 '&NAME'
# =====================================================================

def test_parse_terminator_with_comment():
    """'/ ! end of NEWRUN' 必须结束块, '/' 不得塞进参数数组。"""
    blocks = parse_namelists("&NEWRUN\nH_max=0.001,\n / ! end of NEWRUN\n")
    assert blocks["NEWRUN"] == {"H_max": 0.001}


def test_parse_inline_terminator_with_comment():
    """'H_max=0.001, / ! end' — 行内终止符 + 尾注。"""
    text = "&NEWRUN\nH_max=0.001, / ! end\n&OUTPUT\nZSTOP=1.0,\n/\n"
    blocks = parse_namelists(text)
    assert blocks["NEWRUN"] == {"H_max": 0.001}
    assert blocks["OUTPUT"] == {"ZSTOP": 1.0}


def test_parse_block_name_with_comment():
    """'&NEWRUN ! main block' 块名不得包含注释。"""
    blocks = parse_namelists("&NEWRUN ! main block\nRUN=1,\n/\n")
    assert "NEWRUN" in blocks
    assert blocks["NEWRUN"]["RUN"] == 1


def test_inline_terminator_roundtrip_via_write():
    """带注释终止符的 deck: parse→write→parse 往返不变。"""
    text = "&NEWRUN\nH_max=0.001, / ! end\n/\n"
    blocks = parse_namelists(text)
    assert blocks["NEWRUN"] == {"H_max": 0.001}
    again = parse_namelists(write_namelist("NEWRUN", blocks["NEWRUN"]))
    assert again["NEWRUN"] == {"H_max": 0.001}


def test_format_input_text_inline_terminator():
    """format_input_text 保留行内终止符注释, 输出可再解析。"""
    out = format_input_text("&NEWRUN\nH_max=0.001, / ! end\n/\n")
    blocks = parse_namelists(out)
    assert blocks["NEWRUN"] == {"H_max": 0.001}
    assert "/ ! end" in out


def test_format_namelist_line_terminator_with_comment():
    """format_namelist_line('/ ! end') 保留注释 (直接调用路径)。"""
    assert format_namelist_line("/ ! end of NEWRUN") == " / ! end of NEWRUN"


# =====================================================================
# R2-1-5 (P1-1/P1-2) 元数据补参 + 复数/二维数组写读
# =====================================================================

def test_metadata_dipole_params():
    from astra_tools.deck.metadata import load
    d = load("DIPOLE")
    assert d["n_params"] == 16
    by_name = {p["name"]: p for p in d["params"]}
    for nm in ("D1(,)", "D2(,)", "D3(,)", "D4(,)"):
        p = by_name[nm]
        assert p["type"] == "Double Complex array"
        assert p["unit"] == "m"
        assert p["default"] == "(0.0,0.0)"
    p = by_name["D_Gap(,)"]
    assert p["type"] == "Real*8 array"
    assert p["unit"] == "m"
    assert p["default"] == "0.05"


def test_metadata_modules_module():
    from astra_tools.deck.metadata import load
    d = load("MODULES")
    assert d["n_params"] == 13
    p = next(p for p in d["params"] if p["name"] == "Module(,)")
    assert p["type"] == "Character*40 array"
    assert p["default"] is None


def test_metadata_output_magnetized_params():
    from astra_tools.deck.metadata import load
    d = load("OUTPUT")
    assert d["n_params"] == 35
    by_name = {p["name"]: p for p in d["params"]}
    assert by_name["Lmagnetized"]["default"] == "FALSE"
    assert by_name["Lsub_rot"]["default"] == "FALSE"
    assert by_name["Lmagnetized"]["type"] == "Logical"
    assert by_name["Lsub_rot"]["type"] == "Logical"


def test_metadata_quadrupole_mult_params():
    from astra_tools.deck.metadata import load
    d = load("QUADRUPOLE")
    assert d["n_params"] == 18
    by_name = {p["name"]: p for p in d["params"]}
    assert by_name["Q_mult_a(,)"]["default"] == "0.0"
    assert by_name["Q_mult_b(,)"]["default"] == "0.0"
    assert by_name["Q_mult_a(,)"]["type"] == "Real*8 array"
    assert by_name["Q_xoff()"]["default"] == "0.0"
    assert by_name["Q_yoff()"]["default"] == "0.0"


def test_metadata_charge_merge_params():
    from astra_tools.deck.metadata import load
    d = load("CHARGE")
    assert d["n_params"] == 37
    by_name = {p["name"]: p for p in d["params"]}
    for i in range(1, 11):
        assert "Merge_%d()" % i in by_name, "缺少 Merge_%d" % i
        assert "Integer" in by_name["Merge_%d()" % i]["type"]


def test_metadata_aperture_ap_gr():
    from astra_tools.deck.metadata import load
    d = load("APERTURE")
    assert d["n_params"] == 22
    p = next(p for p in d["params"] if p["name"] == "Ap_GR(,)")
    assert p["type"] == "Real*8 array"
    assert p["unit"] == "mm"
    assert p["default"] is None


def test_write_parse_complex_array_roundtrip():
    """DIPOLE D1=(1.0,2.0) 复数数组写读往返。"""
    text = write_namelist("DIPOLE", {"D1": [1 + 2j, 3 + 4j]})
    assert "D1=(1, 2), (3, 4)," in text
    d = parse_namelists(text)["DIPOLE"]
    assert d["D1"] == [1 + 2j, 3 + 4j]


def test_write_parse_complex_scalar_roundtrip():
    """标量复数写读往返。"""
    text = write_namelist("DIPOLE", {"D1": 1 + 2j})
    d = parse_namelists(text)["DIPOLE"]
    assert d["D1"] == 1 + 2j


def test_write_parse_2d_array_column_major_roundtrip():
    """QUADRUPOLE Q_mult_a 二维数组写读往返 (Fortran 列主序)。"""
    text = write_namelist("QUADRUPOLE", {"Q_mult_a": [[0.1, 0.2], [0.3, 0.4]]})
    assert "Q_mult_a=0.1, 0.3, 0.2, 0.4," in text  # (1,1),(2,1),(1,2),(2,2)
    d = parse_namelists(text)["QUADRUPOLE"]
    assert d["Q_mult_a"] == [0.1, 0.3, 0.2, 0.4]
    # 再写出→再解析, 内容稳定
    text2 = write_namelist("QUADRUPOLE", d)
    assert parse_namelists(text2)["QUADRUPOLE"]["Q_mult_a"] == d["Q_mult_a"]


def test_write_parse_2d_np_array():
    text = write_namelist("DIPOLE", {"D_Gap": np.array([[0.05, 0.06],
                                                        [0.05, 0.06]])})
    d = parse_namelists(text)["DIPOLE"]
    assert d["D_Gap"] == [0.05, 0.05, 0.06, 0.06]


def test_form_getter_complex_array_d1():
    """DIPOLE 表单 D1 文本框 '(1.0, 2.0), (3.0, 4.0)' -> 复数数组。"""
    from astra_tools.widgets.forms import namelist_form
    wmap, getter = namelist_form("DIPOLE", only=["D1"], show=False)
    wmap["D1"].value = "(1.0, 2.0), (3.0, 4.0)"
    out = getter()
    assert out["D1"] == [1 + 2j, 3 + 4j]
    # 写回 deck 可 parse 回读
    text = write_namelist("DIPOLE", out)
    assert parse_namelists(text)["DIPOLE"]["D1"] == [1 + 2j, 3 + 4j]


def test_dipole_metadata_deck_roundtrip(tmp_path):
    """含复数与 2D 参数的整 deck parse→write→parse 往返。"""
    text = ("&DIPOLE\nLDipole=T,\nD1=(1.0, 0.6),\nD2=(-1.0, 0.9),\n"
            "D_Gap=0.05,\n/\n")
    blocks = parse_namelists(text)
    p = write_input_deck(blocks, tmp_path / "d.in")
    assert parse_namelists(p.read_text()) == blocks


# =====================================================================
# R2-1-6 (P2-2) astra.ipynb 表单组参数名 + only= 未匹配警告
# =====================================================================

def _notebook_cell5_source():
    nb = json.loads((PROJECT_ROOT / "notebooks" / "astra.ipynb").read_text())
    cell = nb["cells"][5]
    assert cell["cell_type"] == "code"
    return "".join(cell["source"])


def test_notebook_cell5_parameter_names():
    src = _notebook_cell5_source()
    assert "Wk_filename" in src
    assert "Wk_z" in src
    assert "Ap_R" in src
    for bad in ("File_Wakefield", "W_pos", "AP_radius"):
        assert bad not in src


def test_notebook_group_names_all_exist_in_metadata():
    """cell 5 的 only= 列表全部命中元数据库 (case-insensitive)。"""
    from astra_tools.deck.metadata import load
    src = _notebook_cell5_source()
    groups = {}
    for line in src.splitlines():
        line = line.strip()
        # 只认 'NAME': [ ... ] 形式的表单组行
        if line.startswith('"') and '": [' in line:
            key, _, val = line.partition(":")
            nl = key.strip().strip('"').strip("'")
            names = [n.strip().strip('"').strip("'")
                     for n in val.strip().strip("[],").split(",")]
            if nl and names:
                groups[nl] = names
    assert groups, "cell 5 未解析出表单组"
    for nl, names in groups.items():
        have = {p["name"].rstrip("() ").lower() for p in load(nl)["params"]}
        missing = [n for n in names if n.rstrip("() ").lower() not in have]
        assert not missing, "%s 组参数名不存在: %r" % (nl, missing)


def test_form_only_unknown_param_warns():
    """namelist_form only= 含不存在参数时必须 emit UserWarning。"""
    from astra_tools.widgets.forms import namelist_form
    with pytest.warns(UserWarning, match="NoSuchParam"):
        namelist_form("NEWRUN", only=["NoSuchParam"], show=False)


def test_metadata_param_info_handles_multidim_suffix():
    """param_info 必须能查到 '(,)' 多维数组记法的参数 (R2 收尾).

    旧实现 rstrip("() ") 把 'D1(,)' 归一成 'D1(,' 而查不到。
    """
    from astra_tools.deck.metadata import param_info, _base
    assert _base("D1(,)") == "D1"
    assert _base("MaxE()") == "MaxE"
    assert _base("Ap_GR(,)") == "Ap_GR"
    nl, entry = param_info("D1(,)")
    assert nl == "DIPOLE" and entry["name"] == "D1(,)"
    nl2, entry2 = param_info("Q_mult_a")
    assert nl2 == "QUADRUPOLE"
    nl3, entry3 = param_info("Module")
    assert nl3 == "MODULES"

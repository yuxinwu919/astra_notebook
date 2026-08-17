"""审计修复的回归测试 (纯代码语言层面 bug)。

对应 data/review 阶段之后的对抗性审计报告: 每个测试锁定一个
"退化/边界输入"下此前会静默错值或崩溃、现已修复的行为。
"""

import numpy as np
import pytest


# ---------------------------------------------------------------
# IO 读取器
# ---------------------------------------------------------------

def test_load_columns_single_column_rejected(tmp_path):
    """单列多行 emit 文件应报列数不足, 而非 reshape 成单行出垃圾。"""
    from astra_tools.io.astra_emit import _load_columns
    p = tmp_path / "one_col.txt"
    p.write_text("0.1\n0.2\n0.3\n0.4\n0.5\n0.6\n0.7\n")
    with pytest.raises(ValueError):
        _load_columns(p, 7)


def test_parse_output_file_lowercase_standardizes(tmp_path):
    """小写扩展名 (run.xemit.001) 也应走协方差标准化, 产出 cov_x__xp。"""
    from astra_tools.io.astra_emit import parse_output_file
    p = tmp_path / "run.xemit.001"
    p.write_text("0 0 0 1.0 0 1.0 0\n")
    d = parse_output_file(p)
    assert "cov_x__xp" in d


def test_read_field_single_column_rejected(tmp_path):
    """单列场表应报列数不足, 而非当单点场静默丢弃整张表。"""
    from astra_tools.io.field_map import read_cavity_field, read_solenoid_field
    p = tmp_path / "one_col.dat"
    p.write_text("0.0\n0.1\n0.2\n")
    with pytest.raises(ValueError):
        read_cavity_field(p)
    with pytest.raises(ValueError):
        read_solenoid_field(p)


def test_solenoid_scaled_zero_peak_rejected():
    """全零 Bz 场表 scaled() 应报错, 而非除零崩溃。"""
    from astra_tools.io.field_map import SolenoidField
    s = SolenoidField(z=np.array([0.0, 1.0]), bz0=np.array([0.0, 0.0]))
    with pytest.raises(ValueError):
        s.scaled(0.35)


def test_3d_map_per_value_2x2x2(tmp_path):
    """逐值头网格 (2,2,2) 各维 n=2 时不再被误判成紧凑头 (显式告警)。"""
    from astra_tools.io.field_map import read_3d_field_map
    vals = [2, -0.1, 0.3, 2, -0.2, 0.4, 2, 0.0, 0.5]
    f_true = np.arange(8.0).reshape(2, 2, 2, order="F")
    vals += [float(v) for v in f_true.reshape(-1, order="F")]
    p = tmp_path / "map.dat"
    p.write_text(" ".join(str(v) for v in vals))
    with pytest.warns(UserWarning):
        x, y, z, f = read_3d_field_map(p)
    assert np.array_equal(f, f_true)
    assert np.allclose(x, [-0.1, 0.3]) and np.allclose(z, [0.0, 0.5])


def test_3d_map_compact_2x2x2_ambiguity_documented(tmp_path):
    """(2,2,2) 紧凑头与逐值头 token 完全同形: 固有歧义按逐值头解析并告警。

    紧凑展开本应得 x = [-1.0, 0.0], 但逐值头优先 (手册格式) 会解析成
    [-1.0, 1.0]; 两种解析的数据长度都是 8, 无法从文件内容区分, 因此
    固定按逐值头解析并发出 UserWarning, 不再静默。
    """
    from astra_tools.io.field_map import read_3d_field_map
    vals = [2, -1.0, 1.0, 2, -2.0, 2.0, 2, -3.0, 3.0]
    f_true = np.arange(8.0).reshape(2, 2, 2, order="F")
    vals += [float(v) for v in f_true.reshape(-1, order="F")]
    p = tmp_path / "compact_2.dat"
    p.write_text(" ".join(str(v) for v in vals))
    with pytest.warns(UserWarning):
        x, y, z, f = read_3d_field_map(p)
    assert np.array_equal(f, f_true)
    assert np.allclose(x, [-1.0, 1.0])      # 逐值头: [min, spacing]
    assert np.allclose(y, [-2.0, 2.0])
    assert np.allclose(z, [-3.0, 3.0])


def test_3d_map_truncated_short_file_rejected(tmp_path):
    """4-6 个 token 的截断 3D 场图应报 ValueError, 而非 IndexError。"""
    from astra_tools.io.field_map import read_3d_field_map
    p = tmp_path / "short.dat"
    p.write_text("1 2 3 4 5\n")
    with pytest.raises(ValueError):
        read_3d_field_map(p)


# ---------------------------------------------------------------
# namelist 解析/格式化
# ---------------------------------------------------------------

def test_parse_namelist_array_continuation_line():
    """多行数组的续行 (无 '=' 的行) 应追加, 而非静默丢弃。"""
    from astra_tools.namelist.parse import parse_namelists
    d = parse_namelists("&NEWRUN\nZL=0.0, 1.0,\n  2.0, 3.0\n/\n")
    assert d["NEWRUN"]["ZL"] == [0.0, 1.0, 2.0, 3.0]


def test_parse_namelist_out_of_order_index():
    """乱序索引声明应按索引序 (Fortran 语义), 而非声明序。"""
    from astra_tools.namelist.parse import parse_namelists
    d = parse_namelists("&SOLENOID\nMaxB(2)=20, MaxB(1)=10,\n/\n")
    assert d["SOLENOID"]["MaxB"] == [10, 20]


def test_format_input_text_terminator_with_comment():
    """带尾注的块终止行 '/ ! comment' 应被识别, in_namelist 复位。"""
    from astra_tools.namelist.format import format_input_text
    out = format_input_text("&NEWRUN\nA=1.0\n /  ! end\noutside\n/\n")
    assert "\noutside\n" in out


# ---------------------------------------------------------------
# run/exec + deck
# ---------------------------------------------------------------

def test_run_program_error_title_not_flagged(monkeypatch, tmp_path):
    """exit 0 且 stdout 含 'ERROR STUDY' 标题时应判成功, 不误报。"""
    import subprocess
    import astra_tools.run.exec as exec_mod

    def fake_run(cmd, **kwargs):
        return subprocess.CompletedProcess(cmd, 0, "  ERROR STUDY begins here\n", "")

    monkeypatch.setattr(exec_mod.subprocess, "run", fake_run)
    result = exec_mod.run_program("fake_astra", tmp_path, stream=False)
    assert result.returncode == 0


def test_run_program_real_error_still_flagged(monkeypatch, tmp_path):
    """exit 0 但 stdout 含独立 ERROR 标记时应判失败。"""
    import subprocess
    import pytest
    import astra_tools.run.exec as exec_mod

    def fake_run(cmd, **kwargs):
        return subprocess.CompletedProcess(cmd, 0, "ERROR: bad input\n", "")

    monkeypatch.setattr(exec_mod.subprocess, "run", fake_run)
    with pytest.raises(RuntimeError):
        exec_mod.run_program("fake_astra", tmp_path, stream=False)


def test_solenoid_bz_at_z_nan_field_returns_none(tmp_path):
    """场表含 nan 时应返回 None (降级), 而非静默返回 nan。"""
    from astra_tools.deck.solenoid import solenoid_bz_at_z
    (tmp_path / "sol.dat").write_text("-0.5 1.0\n0.0 nan\n0.5 1.0\n")
    (tmp_path / "astra.in").write_text(
        "&SOLENOID\n LBField=T,\n File_Bfield(1)='sol.dat',\n"
        " MaxB(1)=0.35,\n S_pos(1)=0.0,\n /\n")
    assert solenoid_bz_at_z(tmp_path / "astra.in", 0.0) is None


# ---------------------------------------------------------------
# analysis
# ---------------------------------------------------------------

def test_bff_all_zero_charge_returns_zero():
    """全零电荷束团应返回全零, 而非 0/0 除零产生全 NaN。"""
    from astra_tools.analysis.bff import compute_bff
    r = compute_bff(z=np.array([0.0, 1e-3]), charge=np.array([0.0, 0.0]))
    assert np.all(r.bff == 0.0) and np.all(r.bff_amplitude == 0.0)


# ---------------------------------------------------------------
# plot
# ---------------------------------------------------------------

def test_bff_log_spaced_kmin_zero_rejected():
    """log_spaced k 网格的 kmin 必须 > 0, 否则 log10(0) 使下游崩溃。"""
    from astra_tools.analysis.bff import compute_bff
    with pytest.raises(ValueError):
        compute_bff(z=np.array([0.0, 1e-3]), charge=np.array([1.0, 1.0]),
                    kmin=0.0)


def test_plot_beta_alpha_mismatched_grids():
    """X/Yemit 网格不一致时逐平面插值 bg, 不再广播崩溃。"""
    import matplotlib.pyplot as plt
    from astra_tools.io.astra_emit import EmitData, EmitSet, RefData
    from astra_tools.plot.advanced_plots import plot_beta_alpha

    def _mk(n):
        return EmitData(z=np.linspace(0, 1, n), t=np.zeros(n), avg=np.zeros(n),
                        rms=np.full(n, 1e-3), rmsprime=np.full(n, 1e-3),
                        emit=np.full(n, 1e-6), corr=np.zeros(n))

    emit = EmitSet(x=_mk(80), y=_mk(100), z=_mk(80))
    ref = RefData(z=np.linspace(0, 1, 50), t=np.zeros(50),
                  pz=np.full(50, 1e9), dedz=np.zeros(50), larmor=np.zeros(50),
                  xoff=np.zeros(50), yoff=np.zeros(50), px=np.zeros(50),
                  py=np.zeros(50))
    fig = plot_beta_alpha(emit, ref=ref)
    plt.close(fig)


def test_plot_plasma_profile_zero_column(tmp_path):
    """全零剖面列应降级为空曲线, 而非除零产生 nan。"""
    import matplotlib.pyplot as plt
    from astra_tools.plot.advanced_plots import plot_plasma_profile
    p = tmp_path / "plasma.dat"
    p.write_text("0.0 0.0\n0.1 0.0\n0.2 0.0\n")
    fig = plot_plasma_profile(p, peak_density_cm3=1e18)
    plt.close(fig)


# ---------------------------------------------------------------
# widgets
# ---------------------------------------------------------------

def test_namelist_form_cavity_defaults_getter():
    """全参表单不修改直接 getter() 不应因含括号注释的默认值崩溃。"""
    from astra_tools.widgets.forms import namelist_form
    _, getter = namelist_form("CAVITY", show=False)
    d = getter()
    assert "C_numb" in d


def test_phase_label_no_dot_filename(tmp_path):
    """无点号文件名 (如 README) 应降级返回文件名, 而非 IndexError。"""
    from astra_tools.widgets.selectors import _phase_label
    f = tmp_path / "README"
    f.write_text("hello\n")
    assert _phase_label(f) == "README"


# ---------------------------------------------------------------
# 字体链 (findfont weight 警告)
# ---------------------------------------------------------------

def test_style_font_family_all_have_normal_weight():
    """字体链每个候选字体都应有 normal(400) 字重。

    否则请求 weight=normal 时 findfont 会回退到 300 并发出
    'Failed to find font weight normal, now using 300.' 警告
    (根因: STSong 在 macOS 只有 weight=300 变体, 已从链中移除)。
    """
    from matplotlib import font_manager
    from astra_tools.plot.style import available_font_family
    weights = {}
    for f in font_manager.fontManager.ttflist:
        weights.setdefault(f.name, set()).add(f.weight)
    for name in available_font_family():
        assert 400 in weights.get(name, set()), (
            "字体 %r 缺 normal(400) 字重, 会触发 findfont weight 警告" % name)

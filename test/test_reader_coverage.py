"""Table 4 读取器 L1 单元测试 (批 1a: 补齐 12+ 零覆盖读取函数).

覆盖: read_track_file / read_cathode_file / read_xemit2 / read_tremit /
read_cr_emit / read_larmor / read_density / read_tcheck / read_cemit_file /
parse_field_map_file / expand_tws_field_map / fix_laser_map_header /
format_input_file / backup_directory / get_version / output_file_type。

每条断言列映射 + 单位换算 + 边界 (列数不足/空文件) 错误。
"""
import stat
import sys
from pathlib import Path

import numpy as np
import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.io import astra_misc as misc
from astra_tools.io import field_map as fm
from astra_tools.io.astra_emit import read_cemit_file, output_file_type
from astra_tools.namelist.format import format_input_file
from astra_tools.run.exec import backup_directory, get_version


def _w(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return p


# ---------------- astra_misc 读取器 ----------------

def test_read_track_file_columns_and_units(tmp_path):
    p = _w(tmp_path, "t.track.001",
           "1 5 1.0 2.0 3.0 100.0 200.0 300.0\n")
    d = misc.read_track_file(p)
    assert d["seq"][0] == 1 and d["status"][0] == 5
    assert d["z"][0] == 1.0
    assert d["x"][0] == pytest.approx(2e-3)      # mm -> m
    assert d["y"][0] == pytest.approx(3e-3)
    assert d["Ez"][0] == 100.0 and d["Er"][0] == 200.0 and d["Ey"][0] == 300.0


def test_read_cathode_file_units(tmp_path):
    p = _w(tmp_path, "c.Cathode.001",
           "0.5 2.0 3.0 4.0 5.0 6.0 7.0 8\n")
    d = misc.read_cathode_file(p)
    assert d["z"][0] == 0.5
    assert d["t"][0] == pytest.approx(2e-9)      # ns -> s
    assert d["E_spch"][0] == 3.0 and d["E_acc"][0] == 4.0
    assert d["q"][0] == 5.0 and d["grid_min"][0] == 6.0
    assert d["flag"][0] == 8


def test_read_xemit2_units(tmp_path):
    p = _w(tmp_path, "e.Xemit2.001", "0 1 2 3 4 5 6\n")
    d = misc.read_xemit2(p)
    assert d["K2z"][0] == 1 and d["K3z"][0] == 2
    assert d["eps_red_z"][0] == pytest.approx(3e-6)   # pi mrad mm -> m.rad
    assert d["eps_red_zE"][0] == pytest.approx(6e-6)


def test_read_tremit_units(tmp_path):
    p = _w(tmp_path, "t.TRemit.001", "0 1 2 3 4\n")
    d = misc.read_tremit(p)
    assert d["t"][0] == pytest.approx(1e-9)
    assert d["eps_tr_x"][0] == pytest.approx(2e-6)
    assert d["eps_tr_z"][0] == pytest.approx(4e-6)


def test_read_cr_emit_units(tmp_path):
    p = _w(tmp_path, "c.Cr_emit.001", "0 1 2 3 4 5 6 7\n")
    d = misc.read_cr_emit(p)
    assert d["x_rms"][0] == pytest.approx(2e-3)    # mm -> m
    assert d["eps_x"][0] == pytest.approx(4e-6)
    assert d["q_rest"][0] == 6 and d["q_cross"][0] == 7


def test_read_larmor_units(tmp_path):
    p = _w(tmp_path, "l.Larmor.001", "0 1 2 3\n")
    d = misc.read_larmor(p)
    assert d["t"][0] == pytest.approx(1e-9)
    assert d["avr"][0] == 2 and d["rms"][0] == 3


def test_read_density_units(tmp_path):
    row = "0 1 " + " ".join(str(i) for i in range(2, 12)) + "\n"
    p = _w(tmp_path, "d.Density.001", row)
    d = misc.read_density(p)
    assert d["t"][0] == pytest.approx(1e-9)
    assert d["N"].shape == (1, 5) and d["dens"].shape == (1, 5)
    assert d["N"][0, 0] == 2 and d["dens"][0, 0] == 7


def test_read_tcheck_units(tmp_path):
    p = _w(tmp_path, "t.tcheck.001", "0 1 2 3 4 5 6 7\n")
    d = misc.read_tcheck(p)
    assert d["t"][0] == pytest.approx(1e-9)
    assert d["scaling"].shape == (1, 5)
    assert d["counter"][0] == 7


def test_readers_reject_short_files(tmp_path):
    for name, text in [("t.track.001", "1 2\n"), ("c.Cathode.001", "1\n")]:
        p = _w(tmp_path, name, text)
        with pytest.raises(ValueError):
            misc.read_track_file(p) if "track" in name else misc.read_cathode_file(p)


# ---------------- Cemit 读取器 ----------------

def test_read_cemit_file_columns(tmp_path):
    # 13 列: z + 4x + 4y + 4z; 纵向 keV.mm -> eV.m (x1), 横向 1e-6
    row = "0 " + " ".join(str(v) for v in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]) + "\n"
    p = _w(tmp_path, "c.Cemit.001", row)
    emit = read_cemit_file(p)
    assert emit.emit[0] == pytest.approx(1e-6)        # eps_xn
    assert emit.avg[0] == pytest.approx(2e-6)         # C95 x
    ce = emit._cemit
    assert ce["y"]["eps_n"][0] == pytest.approx(5e-6)
    assert ce["z"]["eps_n"][0] == pytest.approx(9.0)  # keV.mm -> eV.m = x1
    assert ce["z"]["c80"][0] == pytest.approx(12.0)


def test_read_cemit_file_rejects_short(tmp_path):
    p = _w(tmp_path, "c.Cemit.001", "0 1 2\n")
    with pytest.raises(ValueError):
        read_cemit_file(p)


def test_output_file_type():
    assert output_file_type("Example.Xemit.001") == "Xemit"
    assert output_file_type("Example.Yemit.001") == "Yemit"
    assert output_file_type("plain.dat") == "dat"


# ---------------- 场图辅助 ----------------

def test_parse_field_map_1d(tmp_path):
    p = _w(tmp_path, "f.dat", "0 1\n0.1 2\n0.2 3\n")
    attrs, data = fm.parse_field_map_file(p)
    assert attrs["type"] == "astra_1d"
    assert data.shape == (3, 2) and data[-1, 1] == 3


def test_parse_field_map_tws(tmp_path):
    p = _w(tmp_path, "tws.dat", "0.0 0.5 6 3\n0 1\n0.5 2\n")
    attrs, data = fm.parse_field_map_file(p)
    assert attrs["type"] == "astra_tws"
    assert attrs["n"] == 6 and attrs["m"] == 3
    assert data.shape == (2, 2)


def test_expand_tws_field_map():
    z0 = np.linspace(0, 1, 11)
    f0 = np.sin(z0)
    zf, ff = fm.expand_tws_field_map(z0, f0, 0.4, 0.6, 1, 3)
    assert len(zf) == len(ff) and len(zf) > len(z0)
    assert np.all(np.isfinite(ff)) and np.all(np.diff(zf) > 0)


def test_fix_laser_map_header(tmp_path):
    p = _w(tmp_path, "laser.dat",
           "8.1e+01 0.0 1e-3\n9.0e+00 0.0 2e-3\n4.0e+01 0.0 3e-3\n0.1\n")
    out = fm.fix_laser_map_header(p)
    lines = p.read_text().splitlines()
    assert lines[0].split()[0] == "81"
    assert lines[1].split()[0] == "9"      # 中间头行同样取整
    assert lines[2].split()[0] == "40"
    assert lines[3] == "0.1"   # 数据原样保留
    assert len(out) == 3


def test_fix_laser_map_header_rejects_short(tmp_path):
    p = _w(tmp_path, "laser.dat", "1 0 1\n")
    with pytest.raises(ValueError):
        fm.fix_laser_map_header(p)


# ---------------- namelist 格式化 ----------------

def test_format_input_file_normalizes(tmp_path):
    # 语义: True=已归一化/无需处理; False=需要修改 (check_only 不落盘)
    p = _w(tmp_path, "deck.in", "&NEWRUN\r\n Head='x',\r\n RUN=1,\r\n /\r\n")
    assert format_input_file(p, check_only=True) is False
    assert b"\r" in p.read_bytes()         # check_only 不修改 (raw bytes)
    format_input_file(p)
    assert b"\r" not in p.read_bytes()
    assert format_input_file(p, check_only=True) is True   # 已归一化


def test_format_input_file_binary_passthrough(tmp_path):
    # 非 UTF-8 文本 (二进制): 按契约跳过并返回 True, 内容不变
    p = tmp_path / "bin.in"
    p.write_bytes(b"\xff\xfe\x00\x01")
    assert format_input_file(p, check_only=True) is True
    assert p.read_bytes() == b"\xff\xfe\x00\x01"


# ---------------- run.exec 辅助 ----------------

def test_backup_directory(tmp_path):
    src = tmp_path / "src"
    src.mkdir()
    (src / "a.txt").write_text("x")
    root = tmp_path / "backups"
    backup_directory(src, root)
    assert any(root.iterdir())
    assert len(list(root.iterdir())) == 1
    # 第二次备份产生新目录, 不覆盖; 同秒碰撞走 -1 后缀
    backup_directory(src, root)
    dirs = sorted(p.name for p in root.iterdir())
    assert len(dirs) == 2 and dirs[1].endswith("-1")
    assert (root / dirs[1] / "a.txt").read_text() == "x"


def test_get_version(tmp_path):
    script = tmp_path / "tool.sh"
    script.write_text("#!/usr/bin/env python3\nprint('my tool version 9.9.9')\n")
    script.chmod(script.stat().st_mode | stat.S_IEXEC)
    assert get_version(str(script), tmp_path) == "my tool version 9.9.9"


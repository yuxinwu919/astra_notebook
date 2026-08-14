"""后端单元测试 (第 1 层): namelist 读写/格式化, 元数据库, 导出, 可执行文件定位."""

import sys
from pathlib import Path

import numpy as np

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.namelist.write import write_namelist, write_input_deck
from astra_tools.namelist.format import format_input_text
from astra_tools.namelist.parse import parse_namelists
from astra_tools.deck.metadata import NAMELISTS, load, param_info
from astra_tools.run import check_executable


class TestWriteNamelist:
    def test_auto_quote(self):
        t = write_namelist("NEWRUN", {"Head": "test"})
        assert "Head='test'," in t

    def test_already_quoted_kept(self):
        t = write_namelist("NEWRUN", {"Head": "'test'"})
        assert "Head='test'," in t

    def test_bool_tf(self):
        t = write_namelist("N", {"a": True, "b": False})
        assert "a=T," in t and "b=F," in t

    def test_arrays(self):
        t = write_namelist("N", {"v": [1.0, 2.5], "s": ["'a'", "'b'"]})
        assert "v=1, 2.5," in t
        assert "s='a', 'b'," in t

    def test_skip_empty(self):
        t = write_namelist("N", {"a": None, "b": "", "c": [], "d": 1})
        assert "a=" not in t and "b=" not in t and "c=" not in t
        assert "d=1," in t

    def test_input_deck(self, tmp_path):
        p = write_input_deck({"INPUT": {"IPart": 100}, "NEWRUN": {"RUN": 1}},
                             tmp_path / "x.in")
        text = p.read_text()
        assert "&INPUT" in text and "&NEWRUN" in text


class TestFormatInput:
    def test_quotes_preserved(self):
        out = format_input_text("&NEWRUN\n Head = ' a = b '\n /\n")
        assert "' a = b '" in out

    def test_comment_no_comma(self):
        out = format_input_text("&NEWRUN\n RUN=1 ! note\n /\n")
        assert "RUN=1, ! note" in out

    def test_crlf(self):
        out = format_input_text("&N\n a=1,\r\n /\r\n")
        assert "\r" not in out


class TestParseNamelist:
    def test_roundtrip(self):
        d = parse_namelists("&NEWRUN\n RUN=1,\n Track_All=T,\n FNAME='a.ini',\n MaxE(1)=10, MaxE(2)=20,\n /\n")
        assert d["NEWRUN"]["RUN"] == 1
        assert d["NEWRUN"]["Track_All"] is True
        assert d["NEWRUN"]["FNAME"] == "a.ini"
        assert d["NEWRUN"]["MaxE"] == [10, 20]


class TestMetadata:
    def test_all_namelists(self):
        assert len(NAMELISTS) == 14
        for nl in NAMELISTS:
            d = load(nl)
            assert d["n_params"] > 0

    def test_nemit_unit(self):
        nl, p = param_info("Nemit_x")
        assert nl == "INPUT"
        assert "\u03c0" in p["unit"]  # π mm mrad


class TestExport:
    def test_csv_header_units(self, tmp_path):
        from astra_tools.export import export_distribution
        from astra_tools.distribution import Distribution
        d = Distribution.from_arrays(
            np.zeros(5), np.zeros(5), np.zeros(5), np.zeros(5), np.zeros(5),
            np.full(5, 1e6), np.zeros(5), np.full(5, 1e-3),
            status=np.full(5, 5), ref_momentum_eVc=1e6)
        r = export_distribution(d, tmp_path)
        head = "".join(r["csv"].open().readlines()[:3])
        assert "[m]" in head and "[eV/c]" in head and "[nC]" in head
        z = np.load(r["npz"])
        assert len(z["x"]) == 5


class TestExecutable:
    def test_find_astra(self):
        p = check_executable("astra", project_dir=PROJECT_ROOT)
        assert Path(p).exists()

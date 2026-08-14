"""Unit tests for format_input.py."""

import tempfile
from pathlib import Path

import pytest

from format_input import format_input_file, format_file, normalize_line_endings, strip_trailing_spaces


class TestNormalizeLineEndings:
    def test_crlf_to_lf(self):
        assert normalize_line_endings("line1\r\nline2\r\n") == "line1\nline2\n"

    def test_cr_to_lf(self):
        assert normalize_line_endings("line1\rline2") == "line1\nline2"

    def test_lf_unchanged(self):
        assert normalize_line_endings("line1\nline2\n") == "line1\nline2\n"


class TestStripTrailingSpaces:
    def test_trailing_spaces(self):
        assert strip_trailing_spaces("  hello   ") == "  hello"

    def test_trailing_tabs(self):
        assert strip_trailing_spaces("hello\t\t") == "hello"

    def test_no_trailing(self):
        assert strip_trailing_spaces("hello") == "hello"


class TestFormatInputFile:
    def test_simple_namelist(self):
        content = "&NEWRUN\r\n  Head = 'test'\r\n  RUN = 1\r\n /\r\n"
        result = format_input_file(content)
        assert "&NEWRUN" in result
        assert "Head='test'," in result
        assert "RUN=1," in result
        assert "\r" not in result

    def test_preserves_comments(self):
        content = "! This is a comment\n&TEST\n /\n"
        result = format_input_file(content)
        assert "! This is a comment" in result

    def test_preserves_ellipsis(self):
        content = "...\n&TEST\n /\n"
        result = format_input_file(content)
        assert "..." in result

    def test_adds_comma(self):
        content = "&TEST\n  key=value\n /\n"
        result = format_input_file(content)
        assert "key=value," in result

    def test_already_has_comma(self):
        content = "&TEST\n  key=value,\n /\n"
        result = format_input_file(content)
        # Should not double-comma
        assert result.count("key=value,") == 1

    def test_removes_spaces_around_equals(self):
        content = "&TEST\n  key  =  value\n /\n"
        result = format_input_file(content)
        assert "key=value," in result

    def test_normalizes_indentation(self):
        content = "&TEST\n    key=value\n /\n"
        result = format_input_file(content)
        assert "  key=value," in result

    def test_strips_trailing_blank_lines(self):
        content = "&TEST\n /\n\n\n"
        result = format_input_file(content)
        assert not result.endswith("\n\n")

    def test_multiple_namelists(self):
        content = "&FIRST\n  a=1\n /\n&S ECOND\n  b=2\n /\n"
        result = format_input_file(content)
        lines = result.split("\n")
        # Should have blank line between namelists
        assert "&FIRST" in lines
        assert "&S ECOND" in lines


class TestFormatFile:
    def test_in_place_edit(self, tmp_path: Path):
        fpath = tmp_path / "test.in"
        fpath.write_text("&TEST\r\n  key  =  value\r\n /\r\n")
        format_file(fpath)
        result = fpath.read_text()
        assert "key=value," in result
        assert "\r" not in result

    def test_check_only_does_not_modify(self, tmp_path: Path):
        fpath = tmp_path / "test.in"
        original = "&TEST\n  key  =  value\n /\n"
        fpath.write_text(original)
        format_file(fpath, check_only=True)
        assert fpath.read_text() == original

    def test_output_to_different_file(self, tmp_path: Path):
        fpath = tmp_path / "test.in"
        fpath.write_text("&TEST\r\n  key  =  value\r\n /\r\n")
        outpath = tmp_path / "out.in"
        format_file(fpath, output=outpath)
        assert outpath.exists()
        assert "key=value," in outpath.read_text()
        # Original unchanged
        assert "key  =  value" in fpath.read_text()

    def test_already_formatted(self, tmp_path: Path):
        fpath = tmp_path / "test.in"
        content = "&TEST\n  key=value,\n /\n"
        fpath.write_text(content)
        result = format_file(fpath)
        assert result is True  # Already canonical

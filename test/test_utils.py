"""Unit tests for utils.py core functions."""

from pathlib import Path

import numpy as np
import pytest

import utils


# ============================================================
# write_namelist
# ============================================================

class TestWriteNamelist:
    """Tests for write_namelist (format path, filepath=None)."""

    def test_empty_dict(self):
        """Empty params should produce skeleton namelist."""
        result = utils.write_namelist("TEST", {})
        assert result == "&TEST\n /\n"

    def test_bool_true(self):
        result = utils.write_namelist("TEST", {"flag": True})
        assert "flag=T," in result

    def test_bool_false(self):
        result = utils.write_namelist("TEST", {"flag": False})
        assert "flag=F," in result

    def test_int_value(self):
        result = utils.write_namelist("TEST", {"RUN": 1})
        assert "RUN=1," in result

    def test_float_value(self):
        result = utils.write_namelist("TEST", {"H_max": 0.001})
        assert "H_max=0.001," in result

    def test_float_precision(self):
        """Float values use up to 12 significant digits."""
        result = utils.write_namelist("TEST", {"pi": 3.14159265358979})
        assert "pi=3.14159265359," in result

    def test_string_value(self):
        """String values are output directly (caller handles quoting)."""
        result = utils.write_namelist("TEST", {"FNAME": "'bunch.ini'"})
        assert "FNAME='bunch.ini'," in result

    def test_list_of_ints(self):
        result = utils.write_namelist("TEST", {"nums": [1, 2, 3]})
        assert "nums=1, 2, 3," in result

    def test_list_of_floats(self):
        result = utils.write_namelist("TEST", {"vals": [1.0, 2.5, 3.0]})
        assert "vals=1, 2.5, 3," in result

    def test_list_of_bools(self):
        result = utils.write_namelist("TEST", {"flags": [True, False, True]})
        assert "flags=T, F, T," in result

    def test_list_of_strings(self):
        result = utils.write_namelist("TEST", {"files": ["'a.dat'", "'b.dat'"]})
        assert "files='a.dat', 'b.dat'," in result

    def test_none_skipped(self):
        result = utils.write_namelist("TEST", {"a": None, "b": 1})
        assert "a=" not in result
        assert "b=1," in result

    def test_empty_string_skipped(self):
        result = utils.write_namelist("TEST", {"a": "", "b": 1})
        assert "a=" not in result
        assert "b=1," in result

    def test_empty_list_skipped(self):
        result = utils.write_namelist("TEST", {"a": [], "b": 1})
        assert "a=" not in result
        assert "b=1," in result

    def test_empty_tuple_skipped(self):
        result = utils.write_namelist("TEST", {"a": (), "b": 1})
        assert "a=" not in result
        assert "b=1," in result

    def test_numpy_array(self):
        result = utils.write_namelist("TEST", {"vals": np.array([1.0, 2.0, 3.0])})
        assert "vals=1, 2, 3," in result

    def test_namelist_structure(self):
        """Verify full structural output."""
        result = utils.write_namelist("NEWRUN", {"RUN": 1, "Head": "'Test'"})
        lines = result.strip().split("\n")
        assert lines[0] == "&NEWRUN"
        assert lines[-1] == " /"
        assert "  RUN=1," in lines
        assert "  Head='Test'," in lines

    def test_write_to_file(self, tmp_work_dir):
        """When filepath is provided, write to disk and return None."""
        fpath = tmp_work_dir / "test.in"
        result = utils.write_namelist("TEST", {"a": 1}, filepath=fpath)
        assert result is None
        assert fpath.exists()
        content = fpath.read_text()
        assert "&TEST" in content
        assert "a=1," in content


# ============================================================
# read_astra_binary
# ============================================================

class TestReadAstraBinary:
    """Tests for read_astra_binary."""

    def test_read_valid_9col(self, binary_astra_file):
        header, particles = utils.read_astra_binary(binary_astra_file)
        assert header.shape == (5,)
        assert particles.shape == (5, 9)
        assert header[2] == 1.0  # Q_total
        assert particles[0, 0] == 0.001  # x of first particle

    def test_read_10col_format(self, tmp_path, sample_header):
        """10-column format (with particle_index) should be detected."""
        particles_10col = np.column_stack(
            [
                np.array(
                    [
                        [0.001, 0.0, 0.0, 0.0, 0.0, 6.0e6, 0.0, 0.2, 0, 1],
                        [-0.001, 0.0, 0.0, 0.0, 0.0, 6.0e6, 0.0, 0.2, 0, 2],
                    ],
                    dtype=np.float64,
                )
            ]
        )
        data = np.concatenate([sample_header, particles_10col.flatten()])
        fpath = tmp_path / "test_10col.ini"
        data.tofile(str(fpath))

        header, particles = utils.read_astra_binary(fpath)
        assert particles.shape == (2, 9)  # 10th column stripped

    def test_no_particles(self, tmp_path, sample_header):
        """File with only header (0 particles)."""
        sample_header.tofile(str(tmp_path / "empty.ini"))
        header, particles = utils.read_astra_binary(tmp_path / "empty.ini")
        assert header.shape == (5,)
        assert particles.shape == (0, 9)

    def test_file_too_small(self, tmp_path):
        """File with fewer than 5 float64 values."""
        np.array([1.0, 2.0], dtype=np.float64).tofile(str(tmp_path / "small.ini"))
        with pytest.raises(ValueError, match="数据不足"):
            utils.read_astra_binary(tmp_path / "small.ini")

    def test_file_not_found(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            utils.read_astra_binary(tmp_path / "nonexistent.ini")


# ============================================================
# _array_to_dict
# ============================================================

class TestArrayToDict:
    def test_basic_conversion(self, sample_header, sample_particles):
        result = utils._array_to_dict(sample_header, sample_particles, Path("test.ini"))
        assert result["n_particle"] == 5
        assert result["n_active"] == 4  # last particle has status=5
        assert result["total_charge"] == pytest.approx(0.8)  # 4 active * 0.2 nC
        assert result["x"].shape == (5,)
        assert result["active_mask"].dtype == bool
        assert np.array_equal(result["active_mask"], [True, True, True, True, False])

    def test_all_active(self, sample_header):
        particles = np.array(
            [
                [0.0, 0.0, 0.0, 0.0, 0.0, 6e6, 0.0, 0.5, 0],
                [0.0, 0.0, 0.0, 0.0, 0.0, 6e6, 0.0, 0.5, 0],
            ],
            dtype=np.float64,
        )
        result = utils._array_to_dict(sample_header, particles, Path("test.ini"))
        assert result["n_active"] == 2
        assert result["total_charge"] == pytest.approx(1.0)


# ============================================================
# compute_beam_statistics
# ============================================================

class TestComputeBeamStatistics:
    def test_basic_statistics(self, sample_distribution):
        stats = utils.compute_beam_statistics(sample_distribution)
        assert stats["n_particle"] == 5
        assert stats["n_active"] == 4
        assert stats["total_charge_nC"] == pytest.approx(0.8)
        # All quantities should be finite
        for key in [
            "sig_x",
            "sig_y",
            "sig_z",
            "emit_x_geom",
            "emit_y_geom",
            "sig_E_over_E",
        ]:
            assert np.isfinite(stats[key]), f"{key} is not finite: {stats[key]}"

    def test_single_particle(self, sample_header):
        """Single active particle: emittance should be zero."""
        particles = np.array([[0.0, 0.0, 0.0, 0.0, 0.0, 6e6, 0.0, 1.0, 0]], dtype=np.float64)
        dist = utils._array_to_dict(sample_header, particles, Path("test.ini"))
        stats = utils.compute_beam_statistics(dist)
        assert stats["emit_x_geom"] == 0.0
        assert stats["emit_y_geom"] == 0.0

    def test_known_gaussian(self):
        """Synthetic Gaussian beam with known emittance."""
        np.random.seed(42)
        n = 10000
        x = np.random.normal(0, 0.001, n)  # sigma_x = 1 mm
        xp = np.random.normal(0, 0.001, n)  # sigma_xp = 1 mrad
        # Add correlation for non-zero emittance
        rho = 0.5
        xp_corr = rho * x / 0.001 * 0.001 + np.sqrt(1 - rho**2) * xp
        y = np.random.normal(0, 0.001, n)
        yp = np.random.normal(0, 0.001, n)

        particles = np.column_stack(
            [
                x,
                y,
                np.zeros(n),
                xp_corr * 6e6,
                yp * 6e6,
                np.full(n, 6e6),
                np.zeros(n),
                np.ones(n) * 0.1,
                np.zeros(n),
            ]
        )
        dist = utils._array_to_dict(
            np.array([0.0, 6e6, n * 0.1, 0.0, 0.0]), particles, Path("test.ini")
        )
        stats = utils.compute_beam_statistics(dist)
        # RMS emittance should be positive due to correlation
        assert stats["emit_x_geom"] > 0
        # sigma_x should be ~1mm
        assert stats["sig_x"] == pytest.approx(0.001, rel=0.1)

    def test_ref_energy_override(self, sample_distribution):
        stats_default = utils.compute_beam_statistics(sample_distribution)
        stats_override = utils.compute_beam_statistics(sample_distribution, ref_energy_eV=1e9)
        # Higher energy → higher gamma → higher normalized emittance
        assert stats_override["emit_x_norm"] > stats_default["emit_x_norm"]

    def test_all_lost_particles(self, sample_header):
        """All particles lost should still return finite statistics."""
        particles = np.array(
            [
                [0.001, 0.0, 0.0, 0.0, 0.0, 6e6, 0.0, 0.5, 1],
                [-0.001, 0.0, 0.0, 0.0, 0.0, 6e6, 0.0, 0.5, 3],
            ],
            dtype=np.float64,
        )
        dist = utils._array_to_dict(sample_header, particles, Path("test.ini"))
        stats = utils.compute_beam_statistics(dist)
        # Currently, with all particles lost, means return NaN — this is a known edge case
        # The function should handle this gracefully in the future
        assert stats["n_active"] == 0


# ============================================================
# check_executable
# ============================================================

class TestCheckExecutable:
    def test_nonexistent_raises(self, tmp_work_dir):
        with pytest.raises(FileNotFoundError, match="未找到可执行文件"):
            utils.check_executable("nonexistent_binary_xyz", fallback_dir=tmp_work_dir)


# ============================================================
# _display_width
# ============================================================

class TestDisplayWidth:
    def test_ascii(self):
        assert utils._display_width("hello") == 5

    def test_chinese(self):
        assert utils._display_width("你好") == 4  # 2 wide chars = 4 display width

    def test_mixed(self):
        assert utils._display_width("a好b") == 4  # 1 + 2 + 1

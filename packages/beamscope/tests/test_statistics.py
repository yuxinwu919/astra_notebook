"""Tests for beam analysis statistics."""

import numpy as np
import pytest

from beamscope.analysis.statistics import BeamStatistics, compute_statistics, print_statistics


class TestComputeStatistics:
    def test_basic(self, sample_distribution):
        stats = compute_statistics(sample_distribution)
        assert stats.n_particle == 5
        assert stats.n_active == 4
        assert stats.total_charge_nC == pytest.approx(0.8)
        assert stats.sig_x > 0
        assert stats.sig_y > 0
        assert stats.sig_z > 0

    def test_single_particle(self):
        from beamscope.distribution import Distribution
        d = Distribution(
            x=np.array([0.0]),
            y=np.array([0.0]),
            z=np.array([0.0]),
            px=np.array([0.0]),
            py=np.array([0.0]),
            pz=np.array([6e6]),
            clock=np.array([0.0]),
            charge=np.array([1.0]),
            status=np.array([0]),
            ref_momentum_eVc=6e6,
        )
        stats = compute_statistics(d)
        assert stats.emit_x_geom == 0.0
        assert stats.emit_y_geom == 0.0

    def test_all_lost_raises(self):
        from beamscope.distribution import Distribution
        d = Distribution(
            x=np.array([0.0, 0.0]),
            y=np.array([0.0, 0.0]),
            z=np.array([0.0, 0.0]),
            px=np.array([0.0, 0.0]),
            py=np.array([0.0, 0.0]),
            pz=np.array([6e6, 6e6]),
            clock=np.array([0.0, 0.0]),
            charge=np.array([1.0, 1.0]),
            status=np.array([1, 3]),  # all lost
            ref_momentum_eVc=6e6,
        )
        with pytest.raises(ValueError, match="No active particles"):
            compute_statistics(d)

    def test_large_gaussian(self, large_gaussian_distribution):
        stats = compute_statistics(large_gaussian_distribution)
        assert stats.n_active == 1000
        # σ_x should be ~1 mm
        assert stats.sig_x == pytest.approx(0.001, rel=0.15)
        assert stats.sig_y == pytest.approx(0.001, rel=0.15)

    def test_to_dict(self, sample_distribution):
        stats = compute_statistics(sample_distribution)
        d = stats.to_dict()
        assert "sig_x_mm" in d
        assert "emit_x_norm_um" in d
        assert isinstance(d["n_active"], int)

    def test_print_statistics(self, sample_distribution, capsys):
        stats = compute_statistics(sample_distribution)
        print_statistics(stats, title="Test Beam")
        captured = capsys.readouterr()
        assert "Test Beam" in captured.out
        assert "n_active" in captured.out.lower() or "active" in captured.out.lower()


class TestBeamStatistics:
    def test_label(self, sample_distribution):
        stats = compute_statistics(sample_distribution, label="before")
        assert stats.label == "before"

    def test_ref_energy_override(self, sample_distribution):
        stats_low = compute_statistics(sample_distribution, ref_momentum_eVc=6e6)
        stats_high = compute_statistics(sample_distribution, ref_momentum_eVc=1e9)
        # Higher energy → higher gamma → higher normalized emittance
        assert stats_high.emit_x_norm > stats_low.emit_x_norm

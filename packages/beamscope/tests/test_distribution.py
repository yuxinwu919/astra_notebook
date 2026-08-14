"""Tests for Distribution data model."""

import numpy as np
import pytest

from beamscope.distribution import Distribution


class TestDistribution:
    def test_creation(self, sample_distribution):
        d = sample_distribution
        assert d.n_particle == 5
        assert d.n_active == 4
        assert d.is_valid

    def test_active_property(self, sample_distribution):
        d = sample_distribution
        assert np.array_equal(d.active, [True, True, True, True, False])

    def test_active_charge(self, sample_distribution):
        d = sample_distribution
        assert d.active_charge_nC == pytest.approx(0.8)

    def test_filter_active(self, sample_distribution):
        d = sample_distribution
        f = d.filter_active()
        assert f.n_particle == 4
        assert f.n_active == 4
        assert np.all(f.status == 0)

    def test_sample(self, large_gaussian_distribution):
        d = large_gaussian_distribution
        s = d.sample(100, seed=42)
        assert s.n_particle == 100
        assert s.n_active == 100

    def test_from_dict(self, sample_distribution):
        data = {
            "x": sample_distribution.x,
            "y": sample_distribution.y,
            "z": sample_distribution.z,
            "px": sample_distribution.px,
            "py": sample_distribution.py,
            "pz": sample_distribution.pz,
            "clock": sample_distribution.clock / 1e-9,  # s → ns
            "macro_charge": sample_distribution.charge,
            "status_flag": sample_distribution.status.astype(np.float64),
            "header": np.array([0.0, 6e6, 0.9, 0.0, 0.0]),
        }
        d = Distribution.from_dict(data, source="test")
        assert d.n_particle == 5
        # Clock converted from ns → s
        assert d.clock[0] == pytest.approx(0.0)

    def test_repr(self, sample_distribution):
        r = repr(sample_distribution)
        assert "Distribution" in r
        assert "active=4" in r

    def test_summary(self, sample_distribution):
        s = sample_distribution.summary()
        assert "Particles" in s
        assert "MeV" in s

    def test_invalid_distribution(self):
        """Distribution with mismatched array lengths should report is_valid=False."""
        d = Distribution(
            x=np.array([1, 2, 3]),
            y=np.array([1, 2]),  # wrong length
            z=np.array([1, 2, 3]),
            px=np.array([1, 2, 3]),
            py=np.array([1, 2, 3]),
            pz=np.array([1, 2, 3]),
            clock=np.array([1, 2, 3]),
            charge=np.array([1, 2, 3]),
            status=np.array([1, 2, 3]),
        )
        assert not d.is_valid

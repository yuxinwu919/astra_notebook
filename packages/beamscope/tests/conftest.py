"""Shared test fixtures for beamscope."""

from pathlib import Path

import numpy as np
import pytest

from beamscope.distribution import Distribution


@pytest.fixture
def sample_distribution() -> Distribution:
    """A simple 5-particle distribution (4 active, 1 lost)."""
    return Distribution(
        x=np.array([0.001, -0.001, 0.0, 0.002, 0.005]),
        y=np.array([0.0, 0.001, -0.001, 0.002, 0.005]),
        z=np.array([0.0, 0.01, -0.01, 0.02, 0.10]),
        px=np.array([0.0, 100.0, -100.0, 200.0, 0.0]),
        py=np.array([0.0, -50.0, 50.0, 200.0, 0.0]),
        pz=np.array([6.0e6, 6.0e6, 6.0e6, 5.9e6, 6.1e6]),
        clock=np.array([0.0, 1e-10, -1e-10, 2e-10, 1e-9]),
        charge=np.array([0.2, 0.2, 0.2, 0.2, 0.1]),
        status=np.array([0, 0, 0, 0, 5], dtype=np.int32),
        ref_time_ns=0.0,
        ref_momentum_eVc=6.0e6,
        total_charge_nC=0.9,
        source="test.ini",
        format="astra_binary",
    )


@pytest.fixture
def large_gaussian_distribution() -> Distribution:
    """A larger (1000 particle) Gaussian distribution for statistical tests."""
    rng = np.random.default_rng(42)
    n = 1000
    return Distribution(
        x=rng.normal(0, 0.001, n),
        y=rng.normal(0, 0.001, n),
        z=rng.normal(0, 0.002, n),
        px=rng.normal(0, 10, n),
        py=rng.normal(0, 10, n),
        pz=rng.normal(6e6, 1e4, n),
        clock=rng.normal(0, 1e-10, n),
        charge=np.ones(n) * 0.1,
        status=np.zeros(n, dtype=np.int32),
        ref_momentum_eVc=6e6,
        total_charge_nC=n * 0.1,
    )


@pytest.fixture
def binary_astra_file(tmp_path: Path, sample_distribution: Distribution) -> Path:
    """Create a synthetic ASTRA binary file for testing I/O."""
    header = np.array([0.0, 6e6, 0.9, 0.0, 0.0], dtype=np.float64)
    particles = np.column_stack([
        sample_distribution.x,
        sample_distribution.y,
        sample_distribution.z,
        sample_distribution.px,
        sample_distribution.py,
        sample_distribution.pz,
        sample_distribution.clock / 1e-9,  # s → ns for ASTRA format
        sample_distribution.charge,
        sample_distribution.status.astype(np.float64),
    ])
    data = np.concatenate([header, particles.flatten()])
    fpath = tmp_path / "test_bunch.ini"
    data.tofile(str(fpath))
    return fpath

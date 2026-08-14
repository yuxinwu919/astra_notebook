"""Shared pytest fixtures for ASTRA Notebook tests."""

import sys
from pathlib import Path

import numpy as np
import pytest

# Add project root to sys.path so tests can import utils, format_input, etc.
PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


@pytest.fixture(scope="session")
def project_root() -> Path:
    """Project root directory."""
    return PROJECT_ROOT


@pytest.fixture
def tmp_work_dir(tmp_path: Path) -> Path:
    """Temporary working directory for I/O tests."""
    return tmp_path


@pytest.fixture
def sample_header() -> np.ndarray:
    """Synthetic ASTRA binary header: ref_time [ns], ref_energy [eV], Q_total [nC],
    reserved, reserved."""
    return np.array([0.0, 6.0e6, 1.0, 0.0, 0.0], dtype=np.float64)


@pytest.fixture
def sample_particles() -> np.ndarray:
    """Synthetic 9-column particle data (5 particles).

    Columns: x[m], y[m], z[m], px[eV/c], py[eV/c], pz[eV/c],
             clock[ns], macro_charge[nC], status_flag
    """
    return np.array(
        [
            # x, y, z, px, py, pz, clock, charge, status
            [0.001, 0.000, 0.00, 0.0, 0.0, 6.0e6, 0.0, 0.2, 0],
            [-0.001, 0.001, 0.01, 100.0, -50.0, 6.0e6, 0.1, 0.2, 0],
            [0.000, -0.001, -0.01, -100.0, 50.0, 6.0e6, -0.1, 0.2, 0],
            [0.002, 0.002, 0.02, 200.0, 200.0, 5.9e6, 0.2, 0.2, 0],
            [0.005, 0.005, 0.10, 0.0, 0.0, 6.1e6, 1.0, 0.1, 5],  # lost particle
        ],
        dtype=np.float64,
    )


@pytest.fixture
def sample_distribution(sample_header: np.ndarray, sample_particles: np.ndarray) -> dict:
    """Build a standardized distribution dict (matching utils._array_to_dict output)."""
    header = sample_header
    x, y, z = sample_particles[:, 0], sample_particles[:, 1], sample_particles[:, 2]
    px, py, pz = sample_particles[:, 3], sample_particles[:, 4], sample_particles[:, 5]
    clock = sample_particles[:, 6]
    charge = sample_particles[:, 7]
    status = sample_particles[:, 8]
    active = status == 0

    return {
        "filepath": Path("synthetic.ini"),
        "header": header,
        "x": x,
        "y": y,
        "z": z,
        "px": px,
        "py": py,
        "pz": pz,
        "clock": clock,
        "macro_charge": charge,
        "status_flag": status,
        "n_particle": len(sample_particles),
        "n_active": int(np.sum(active)),
        "total_charge": float(np.sum(charge[active])),
        "active_mask": active,
    }


@pytest.fixture
def binary_astra_file(tmp_path: Path, sample_header: np.ndarray, sample_particles: np.ndarray) -> Path:
    """Write a synthetic ASTRA binary file to a temp directory."""
    data = np.concatenate([sample_header, sample_particles.flatten()])
    filepath = tmp_path / "test_bunch.ini"
    data.tofile(str(filepath))
    return filepath

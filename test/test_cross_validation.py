"""Golden-sample cross-validation against real ASTRA output.

Validates the entire statistics chain (reader + physics) against the
official ASTRA Manual Example output files:

    examples/Manual_Example/Example.0150.001   (ASCII phase space, z=1.5 m)
    examples/Manual_Example/Example.Xemit.001  (transverse emit evolution)
    examples/Manual_Example/Example.Zemit.001  (longitudinal emit evolution)
    examples/Manual_Example/Solenoid.dat       (on-axis field table)

Reference values: last row of Xemit/Zemit. Agreement expected < 0.5%.

The example beamline contains a solenoid (MaxB=0.35 T at S_pos=1.2 m),
so the emittance MUST be computed with the canonical momentum
(p~x = px + c*Bz*y/2, manual 4.13.1).
"""

from pathlib import Path

import numpy as np
import pytest
from scipy.interpolate import interp1d

import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from astra_tools.io import read_distribution
from astra_tools.analysis.statistics import compute_statistics

DATA = PROJECT_ROOT / "examples" / "Manual_Example"


def _bz_at_bunch_center() -> float:
    """On-axis solenoid field at the bunch position z=1.5 m."""
    table = np.loadtxt(DATA / "Solenoid.dat")
    maxB = 0.35          # Example.in: MaxB(1)=0.35
    s_pos = 1.2          # Example.in: S_pos(1)=1.2
    z_bunch = 1.5        # ZSTOP
    value = interp1d(table[:, 0], table[:, 1])(z_bunch - s_pos)
    return float(value * maxB / table[:, 1].max())


@pytest.fixture(scope="module")
def stats():
    dist = read_distribution(DATA / "Example.0150.001")
    return compute_statistics(dist, bz_on_axis_T=_bz_at_bunch_center())


def test_reader_header(stats):
    """Reference particle quantities recovered from the first row."""
    assert stats.ref_momentum_eVc == pytest.approx(1.0005e9, rel=1e-4)
    assert stats.ref_kinetic_energy_eV == pytest.approx(999.99e6, rel=1e-4)


def test_transverse_moments_match_xemit(stats):
    xemit = np.loadtxt(DATA / "Example.Xemit.001")
    last = xemit[-1]
    assert stats.mean_x * 1e3 == pytest.approx(last[2], rel=5e-3)   # mm
    assert stats.sig_x * 1e3 == pytest.approx(last[3], rel=5e-4)    # mm
    assert stats.sig_xp * 1e3 == pytest.approx(last[4], rel=5e-4)   # mrad
    # Xemit column 6 stores eps_n in units of 1e-6 m.rad (mm mrad)
    assert stats.emit_x_norm * 1e6 == pytest.approx(last[5], rel=5e-3)


def test_longitudinal_moments_match_zemit(stats):
    zemit = np.loadtxt(DATA / "Example.Zemit.001")
    last = zemit[-1]
    assert stats.mean_E_kin_eV * 1e-6 == pytest.approx(last[2], rel=1e-4)  # MeV
    assert stats.sig_z * 1e3 == pytest.approx(last[3], rel=1e-3)           # mm
    assert stats.sig_E_eV * 1e-3 == pytest.approx(last[4], rel=1e-3)       # keV
    # longitudinal emittance (uncorrelated): sigma_E * sigma_z
    eps_zn = stats.sig_E_eV * 1e-3 * stats.sig_z * 1e3
    assert eps_zn == pytest.approx(last[5], rel=2e-3)


def test_reference_particle_is_absolute():
    """The reference row (z=1.5) must not pollute the bunch statistics."""
    dist = read_distribution(DATA / "Example.0150.001")
    # Bunch z values are absolute around 1.5 m, not 0
    assert np.mean(dist.z) == pytest.approx(1.5, abs=0.01)
    assert dist.ref_z_m == pytest.approx(1.5)

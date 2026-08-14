"""ASTRA field-map readers: cavity E-fields and solenoid fields.

File formats (ASTRA Manual V3.2, sections 6.9 / 6.10):

  * Cavity E-field table (File_Efield): two columns, free format:
        z [m], Ez [MV/m]        (longitudinal on-axis field)
    ASTRA scales the table to MaxE(n) unless C_noscale=T; the raw file
    values are in MV/m. Transverse and azimuthal components follow from
    the axis field by Maxwell's equations (off-axis expansion).

  * Solenoid field table (File_Bfield): two columns, free format:
        z [m], Bz [arb. units]  (longitudinal on-axis field)
    Scaled to MaxB(n) unless S_noscale=T.

  * Wake potential table (File_Wakefield, section 6.8):
        first line: N 0, then N lines: s [m], W [V/C]

Off-axis field expansion (used by fieldplot, manual chapter 8):
    Ez(r,z) = Ez0 - (r^2/4) Ez0'' + (r^4/64) Ez0''''
    Er(r,z) = -(r/2) Ez0' + (r^3/16) Ez0'''
    Bphi(r,z) = (r/(2 c^2)) dEz0/dt      (RF: dEz0/dt = omega * Ez0)
    Br(r,z)  = -(r/2) Bz0'
    Bz(r,z)  = Bz0 - (r^2/4) Bz0''
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..constants import C_LIGHT


@dataclass
class CavityField:
    """On-axis cavity field table.

    z [m], ez0 [V/m] (raw file values in MV/m converted to V/m).
    """

    z: np.ndarray
    ez0: np.ndarray
    source: str = ""

    @property
    def peak_field_MVpm(self) -> float:
        return float(np.max(np.abs(self.ez0))) * 1e-6

    def field_at(self, r: np.ndarray, z: np.ndarray, omega: float = 0.0):
        """Off-axis field components at (r, z) via the axis expansion.

        Args:
            r: radial positions [m] (array).
            z: longitudinal positions [m] (array, same shape as r).
            omega: RF angular frequency [rad/s]; 0 for static fields
                   (Bphi = 0).

        Returns:
            (ez, er, bphi) arrays [V/m], [V/m], [T].
        """
        r = np.asarray(r, dtype=float)
        z = np.asarray(z, dtype=float)
        ez0 = np.interp(z, self.z, self.ez0)
        dez = np.gradient(self.ez0, self.z)
        d2ez = np.gradient(dez, self.z)
        d3ez = np.gradient(d2ez, self.z)
        d4ez = np.gradient(d3ez, self.z)
        ez0p = np.interp(z, self.z, dez)
        ez0pp = np.interp(z, self.z, d2ez)
        ez0ppp = np.interp(z, self.z, d3ez)
        ez0pppp = np.interp(z, self.z, d4ez)

        ez = ez0 - (r**2 / 4.0) * ez0pp + (r**4 / 64.0) * ez0pppp
        er = -(r / 2.0) * ez0p + (r**3 / 16.0) * ez0ppp
        bphi = (r / (2.0 * C_LIGHT**2)) * (omega * ez0) if omega else np.zeros_like(r)
        return ez, er, bphi


@dataclass
class SolenoidField:
    """On-axis solenoid field table.

    z [m], bz0 [arb. units]; scale to physical Tesla with the MaxB value.
    """

    z: np.ndarray
    bz0: np.ndarray
    source: str = ""

    def scaled(self, max_b_T: float) -> "SolenoidField":
        """Return a copy scaled so that peak |Bz| equals max_b_T [T]."""
        peak = float(np.max(np.abs(self.bz0)))
        out = SolenoidField(z=self.z.copy(), bz0=self.bz0 * (max_b_T / peak),
                            source=self.source)
        return out

    def field_at(self, r: np.ndarray, z: np.ndarray):
        """Off-axis (Br, Bz) via the axis expansion.

        Args:
            r: radial positions [m].
            z: longitudinal positions [m].

        Returns:
            (br, bz) [T].
        """
        r = np.asarray(r, dtype=float)
        z = np.asarray(z, dtype=float)
        bz0 = np.interp(z, self.z, self.bz0)
        dbz = np.gradient(self.bz0, self.z)
        d2bz = np.gradient(dbz, self.z)
        bz0p = np.interp(z, self.z, dbz)
        bz0pp = np.interp(z, self.z, d2bz)
        br = -(r / 2.0) * bz0p
        bz = bz0 - (r**2 / 4.0) * bz0pp
        return br, bz


@dataclass
class WakePotential:
    """Wake potential table: s [m], W [V/C] (pseudo-Green function)."""

    s: np.ndarray
    w: np.ndarray
    source: str = ""


def read_cavity_field(path) -> CavityField:
    """Read a two-column cavity field table (z [m], Ez [MV/m])."""
    path = Path(path)
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError("cavity field file must have >= 2 columns: " + str(path))
    return CavityField(z=data[:, 0].astype(float),
                       ez0=data[:, 1].astype(float) * 1e6,  # MV/m -> V/m
                       source=str(path))


def read_solenoid_field(path) -> SolenoidField:
    """Read a two-column solenoid field table (z [m], Bz [arb. units])."""
    path = Path(path)
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError("solenoid field file must have >= 2 columns: " + str(path))
    return SolenoidField(z=data[:, 0].astype(float),
                         bz0=data[:, 1].astype(float),
                         source=str(path))


def read_wake_potential(path) -> WakePotential:
    """Read an ASTRA wake table: header line 'N 0' then N rows (s, W).

    Multiple blocks (e.g. monopole + dipole) may follow each other;
    only the first block is returned. See manual section 6.8.
    """
    path = Path(path)
    lines = [ln.split() for ln in path.read_text().splitlines()
             if ln.strip() and not ln.strip().startswith(("#", "!"))]
    if not lines:
        raise ValueError("empty wake file: " + str(path))
    n = int(float(lines[0][0]))
    rows = np.asarray([[float(v) for v in ln] for ln in lines[1:1 + n]], dtype=float)
    return WakePotential(s=rows[:, 0], w=rows[:, 1], source=str(path))

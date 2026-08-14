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


# ============================================================
# 1D field maps incl. travelling-wave structures (TWS)
# ============================================================
# Vendored and adapted from lume-astra (C. Mayes et al., LBNL,
# Apache-2.0): https://github.com/ChristopherMayes/lume-astra
# (astra/fieldmaps.py). ASTRA manual V3.2 section 6.9 describes the
# TWS field-map format: a 4-value header line 'z1 z2 n m' followed by
# a two-column (z, Ez) body.


def parse_field_map_file(path):
    """Parse a 1D ASTRA field map (standard 2-column or TWS).

    Returns (attrs, data):
        attrs = {'type': 'astra_1d' | 'astra_tws', 'z1','z2','n','m'}
        data  = (N, 2) array [z, field]
    """
    path = Path(path)
    with open(path) as f:
        header = [float(v) for v in f.readline().split()]

    attrs = {}
    if len(header) == 4:
        attrs["type"] = "astra_tws"
        attrs["z1"] = header[0]
        attrs["z2"] = header[1]
        attrs["n"] = int(header[2])
        attrs["m"] = int(header[3])
        data = np.loadtxt(path, skiprows=1)
    else:
        attrs["type"] = "astra_1d"
        data = np.loadtxt(path)
    return attrs, data


def expand_tws_field_map(z0, f0, z1, z2, m_cells_in_body, n_cell):
    """Periodically expand a TWS field-map body over n_cell cells.

    Layout: |Entrance| Cells | Exit | -> |Entrance| Cells |...|Cells| Exit |
    (vendored from lume-astra's expand_tws_fmap).

    Args:
        z0: z positions of the raw table [m].
        f0: field values of the raw table.
        z1, z2: start/end of the periodic cell body [m].
        m_cells_in_body: number of cells inside the raw body (attrs 'm').
        n_cell: total number of cells requested (C_numb).

    Returns:
        (zfull, ffull) arrays.
    """
    z0 = np.asarray(z0, dtype=float)
    f0 = np.asarray(f0, dtype=float)
    zmin, zmax = float(z0.min()), float(z0.max())
    dz = float(np.mean(np.diff(z0)))
    l_entrance = z1 - zmin
    l_exit = zmax - z2
    l_cell = z2 - z1

    n_repeat = int(n_cell / m_cells_in_body)

    z_entrance = np.linspace(zmin, z1, int(round(l_entrance / dz + 1)))
    z_cell = np.linspace(z1, z2, int(round(l_cell / dz + 1)))
    z_exit = np.linspace(z2, zmax, int(round(l_exit / dz + 1)))

    f_entrance = np.interp(z_entrance, z0, f0)
    f_cell = np.interp(z_cell, z0, f0)
    f_exit = np.interp(z_exit, z0, f0)

    ztot = [z_entrance[:-1]]
    ftot = [f_entrance[:-1]]
    for i in range(n_repeat):
        ztot.append(z_cell[:-1] + i * l_cell)
        ftot.append(f_cell[:-1])
    ztot.append(z_exit + (n_repeat - 1) * l_cell)
    ftot.append(f_exit)

    return np.concatenate(ztot), np.concatenate(ftot)


def fix_laser_map_header(path):
    """Convert a MATLAB-written laser 3D map to ASTRA grid-header format.

    DESY ships Plasma_Example_2 laser.dat with three header lines in
    (n, min, spacing) form; the ASTRA 3D map reader requires explicit
    grid value lines (n followed by the n values). Converts in place.
    """
    path = Path(path)
    with open(path) as f:
        lines = f.readlines()
    if len(lines) < 4:
        raise ValueError("laser map too short: " + str(path))
    grids = []
    for i in range(3):
        vals = [float(v) for v in lines[i].split()]
        n = int(vals[0])
        x0, dx = vals[1], vals[2]
        grid = [x0 + j * dx for j in range(n)]
        grids.append((n, grid))
    with open(path, "w") as f:
        for n, grid in grids:
            f.write(" ".join([str(n)] + ["%.16e" % v for v in grid]) + "\n")
        f.writelines(lines[3:])
    return [(n, g[0], g[-1]) for n, g in grids]


def read_3d_field_map(path):
    """读取 ASTRA 3D 场图 (如 3D_test.ex / 3D_Dipole.bx).

    格式: 3 个网格行 (n 与 n 个网格值, 自由格式), 随后为数据值,
    顺序为 x 最快、y 次之、z 最慢 (Fortran 序)。

    Returns:
        (x, y, z, F) — F 为 (nx, ny, nz) 数组, 按 F[ix, iy, iz] 索引,
        SI 单位按文件名约定 (ex/ey/ez: V/m; bx/by/bz: T; 数值原样返回)。
    """
    path = Path(path)
    toks = path.read_text().split()
    vals = [float(t) for t in toks]
    idx = 0

    def read_grid():
        nonlocal idx
        n = int(vals[idx])
        idx += 1
        g = np.array(vals[idx:idx + n], dtype=float)
        idx += n
        return g

    x = read_grid()
    y = read_grid()
    z = read_grid()
    nx, ny, nz = len(x), len(y), len(z)
    data = np.array(vals[idx:idx + nx * ny * nz], dtype=float)
    if len(data) < nx * ny * nz:
        raise ValueError("3D map data truncated: " + str(path))
    f = data.reshape(nx, ny, nz, order="F")  # x fastest -> F[ix, iy, iz]
    return x, y, z, f

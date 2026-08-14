"""ASTRA evolution-file readers: Xemit/Yemit/Zemit/Cemit/Sigma/ref/Log.

All formats per ASTRA Manual V3.2, Table 4. Units are converted to SI
on read; the original display units are documented per column.

    Xemit / Yemit   (7 cols): z[m] t[ns] u_avr[mm] u_rms[mm]
                              u'_rms[mrad] eps_n[1e-6 m.rad] <uu'>_avr[mm.mrad]
    Zemit           (7 cols): z[m] t[ns] E_kin[MeV] z_rms[mm]
                              dE_rms[keV] eps_zn[keV.mm] <z E'>_avr[keV.mm]
    ref             (9 cols): z[m] t[ns] pz[MeV/c] dE/dz[MeV/m]
                              Larmor[rad] x_off[mm] y_off[mm] px[eV/c] py[eV/c]
    Sigma           (23 cols): z[m] E_kin[MeV] + 21 upper-triangle elements
                              of the 6x6 covariance matrix in canonical
                              coordinates (x, p~x, y, p~y, z, E_kin)
    Cemit           (13 cols): z[m] + 12 core emittances
                              eps_xn, Cx95, Cx90, Cx80 [1e-6 m.rad],
                              eps_yn, Cy95, Cy90, Cy80 [1e-6 m.rad],
                              eps_zn, Cz95, Cz90, Cz80 [keV.mm]

Emittance unit convention (validated against ASTRA output, lume-astra
parsers and pmd-beamphysics):
  * ASTRA prints emittances in 'pi mm mrad' / 'pi keV mm'. The pi marks
    the quantity as the AREA of the RMS phase-space ellipse
    (area = pi*a*b); the RMS statistical emittance from particle data is
    eps_rms = a*b. Numerically 'pi mm mrad' and 'mm mrad' are therefore
    the SAME value in ASTRA files: file value x 1e-6 -> m.rad, do NOT
    multiply by pi. lume-astra labels the column 'mm-mrad' (factor 1e-6)
    and pmd-beamphysics computes norm_emit = sqrt(det cov(x,px))/mc2
    in meters - all consistent.
  * the last column stores cov(u, u') / sigma_u (lume-astra name:
    cov_x__xp/sigma_x), in mrad.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)

MM_TO_M = 1e-3
MRAD_TO_RAD = 1e-3
MEV_TO_EV = 1e6
KEV_TO_EV = 1e3
NS_TO_S = 1e-9


def _load_columns(path: Path, ncols: int) -> np.ndarray:
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < ncols:
        raise ValueError(
            "%s: expected >= %d columns, got %d" % (path.name, ncols, data.shape[1])
        )
    return data


def _base_from_rootname(rootname: str) -> str:
    root = Path(rootname)
    return str(root.parent / root.stem) if root.suffix else str(root)


# ============================================================
# Data models
# ============================================================

@dataclass
class EmitData:
    """Single-plane emittance evolution (7 columns, SI units).

    For Xemit/Yemit: avg=<u> [m], rms=sigma_u [m], rmsprime=sigma_u' [rad],
    emit=eps_n [m.rad], corr=<u u'> [m.rad].
    For Zemit: avg=E_kin [eV], rms=sigma_z [m], rmsprime=sigma_E [eV],
    emit=eps_zn [eV.m], corr=<z E'> [eV.m].
    """

    z: np.ndarray          # [m]
    t: np.ndarray          # [s]
    avg: np.ndarray
    rms: np.ndarray
    rmsprime: np.ndarray
    emit: np.ndarray
    corr: np.ndarray
    label: str = ""        # 'x', 'y', 'z'
    plane: str = ""        # 'horizontal', 'vertical', 'longitudinal'


@dataclass
class EmitSet:
    """Full set of X/Y/Z emittance evolution + optional Cemit."""

    x: EmitData
    y: EmitData
    z: EmitData
    c: Optional[EmitData] = None
    filename: str = ""


@dataclass
class SigmaData:
    """6x6 beam covariance matrix evolution (23 columns).

    Covariance in canonical coordinates (x, p~x, y, p~y, z, E_kin),
    manual 4.13. matrix[i,j] in SI: positions [m^2] on (1,1),(3,3),(5,5),
    mixed terms [m.eV/c] etc.
    """

    z: np.ndarray          # [m]
    e_kin_eV: np.ndarray   # [eV] (file stores MeV)
    matrix: np.ndarray     # (N, 6, 6) covariance
    enx: np.ndarray        # eigen-emittance x [m.rad]
    eny: np.ndarray        # eigen-emittance y [m.rad]
    enz: np.ndarray        # longitudinal emittance [eV.m]
    filename: str = ""


@dataclass
class RefData:
    """Reference particle trajectory (9 columns, SI units)."""

    z: np.ndarray       # [m]
    t: np.ndarray       # [s]
    pz: np.ndarray      # [eV/c]
    dedz: np.ndarray    # [eV/m]
    larmor: np.ndarray  # [rad]
    xoff: np.ndarray    # [m]
    yoff: np.ndarray    # [m]
    px: np.ndarray      # [eV/c]
    py: np.ndarray      # [eV/c]
    filename: str = ""


# ============================================================
# Readers
# ============================================================

def read_emit_files(rootname: str, run: str = "001") -> EmitSet:
    """Read Xemit/Yemit/Zemit (and optionally Cemit) for one run."""
    base = _base_from_rootname(rootname)
    x = _read_emit_plane(Path(base + ".Xemit." + run), "x", "horizontal",
                         u_is_energy=False)
    y = _read_emit_plane(Path(base + ".Yemit." + run), "y", "vertical",
                         u_is_energy=False)
    z = _read_emit_plane(Path(base + ".Zemit." + run), "z", "longitudinal",
                         u_is_energy=True)

    c = None
    cpath = Path(base + ".Cemit." + run)
    if cpath.exists():
        try:
            c = read_cemit_file(cpath)
        except Exception as e:
            logger.warning("could not read Cemit file: %s", e)

    name = Path(base).name
    logger.info(
        "Read emit files '%s': Xemit=%d, Yemit=%d, Zemit=%d rows",
        name, len(x.z), len(y.z), len(z.z),
    )
    return EmitSet(x=x, y=y, z=z, c=c, filename=name)


def _read_emit_plane(path: Path, label: str, plane: str, u_is_energy: bool) -> EmitData:
    if not path.exists():
        raise FileNotFoundError("emit file not found: " + str(path))
    data = _load_columns(path, 7)
    z = data[:, 0]
    t = data[:, 1] * NS_TO_S

    if u_is_energy:
        avg = data[:, 2] * MEV_TO_EV        # E_kin [eV]
        rms = data[:, 3] * MM_TO_M          # z_rms [m]
        rmsprime = data[:, 4] * KEV_TO_EV   # dE_rms [eV]
        emit = data[:, 5] * (KEV_TO_EV * MM_TO_M)  # eps_zn [eV.m]
        corr = data[:, 6] * (KEV_TO_EV * MM_TO_M)  # cov(z,E)/sigma_z [eV.m]
    else:
        avg = data[:, 2] * MM_TO_M          # <u> [m]
        rms = data[:, 3] * MM_TO_M          # sigma_u [m]
        rmsprime = data[:, 4] * MRAD_TO_RAD  # sigma_u' [rad]
        emit = data[:, 5] * 1e-6            # eps_n [m.rad] (no pi factor!)
        corr = data[:, 6] * MRAD_TO_RAD     # cov(u,u')/sigma_u [rad]

    return EmitData(z=z, t=t, avg=avg, rms=rms, rmsprime=rmsprime,
                    emit=emit, corr=corr, label=label, plane=plane)


def read_cemit_file(path) -> EmitData:
    """Read a Cemit (core emittance) file: 13 columns.

    z + (eps_xn, Cx95, Cx90, Cx80) + (eps_yn, Cy95, Cy90, Cy80)
      + (eps_zn, Cz95, Cz90, Cz80)
    Transverse core emittances in 1e-6 m.rad units, longitudinal in
    keV.mm. Stored in the generic EmitData model: avg/rms/rmsprime carry
    C95/C90/C80 (x-plane), corr unused.
    """
    path = Path(path)
    data = _load_columns(path, 13)
    z = data[:, 0]
    # x plane: eps_n at cols 1..4; store C95/C90/C80 in avg/rms/rmsprime
    emit_x = data[:, 1] * 1e-6   # eps_xn [m.rad]
    c95_x = data[:, 2] * 1e-6
    c90_x = data[:, 3] * 1e-6
    c80_x = data[:, 4] * 1e-6
    # y plane
    emit_y = data[:, 5] * 1e-6
    c95_y = data[:, 6] * 1e-6
    c90_y = data[:, 7] * 1e-6
    c80_y = data[:, 8] * 1e-6
    # z plane (keV.mm)
    emit_z = data[:, 9] * (KEV_TO_EV * MM_TO_M)
    c95_z = data[:, 10] * (KEV_TO_EV * MM_TO_M)
    c90_z = data[:, 11] * (KEV_TO_EV * MM_TO_M)
    c80_z = data[:, 12] * (KEV_TO_EV * MM_TO_M)

    emit = EmitData(
        z=z, t=np.zeros(len(z)),
        avg=c95_x, rms=c90_x, rmsprime=c80_x,
        emit=emit_x, corr=np.zeros(len(z)),
        label="core_x", plane="horizontal",
    )
    # Store the other planes in attrs-like arrays via a dict attribute
    emit._cemit = {
        "x": {"eps_n": emit_x, "c95": c95_x, "c90": c90_x, "c80": c80_x},
        "y": {"eps_n": emit_y, "c95": c95_y, "c90": c90_y, "c80": c80_y},
        "z": {"eps_n": emit_z, "c95": c95_z, "c90": c90_z, "c80": c80_z},
    }
    return emit


def read_sigma_file(rootname: str, run: str = "001") -> SigmaData:
    """Read a Sigma (6x6 covariance) file: 23 columns.

    z[m], E_kin[MeV], then the 21 upper-triangle elements sig(i,j),
    i=1..6, j=i..6 (row-major) of the covariance matrix in canonical
    coordinates (x, p~x, y, p~y, z, E_kin) - manual 4.13.

    Empirical validation (examples/Manual_Example, last row):
      * sig(1,1) [m^2] = sigma_x^2 and sig(5,5) [m^2] = sigma_z^2
        (exact match with Xemit/Zemit)
      * sig(2,2) and sig(6,6) scale like sigma_x'^2 and sigma_E^2 but
        with an additional factor of ~3.83 whose origin is not documented
        in the manual; the eigen-emittances derived here are therefore
        marked experimental. For validated emittances use the Xemit
        files (read_emit_files) instead.
    """
    base = _base_from_rootname(rootname)
    path = Path(base + ".Sigma." + run)
    data = _load_columns(path, 23)

    z = data[:, 0]
    e_kin = data[:, 1] * MEV_TO_EV
    cov_elements = data[:, 2:]  # (N, 21)

    n = len(z)
    matrix = np.zeros((n, 6, 6))
    idx = 0
    for i in range(6):
        for j in range(i, 6):
            matrix[:, i, j] = cov_elements[:, idx]
            if i != j:
                matrix[:, j, i] = cov_elements[:, idx]
            idx += 1

    # Eigen-emittances: imaginary parts of eigenvalues of Sigma @ J
    # (4x4 transverse block), plus longitudinal sqrt(det).
    enx = np.zeros(n)
    eny = np.zeros(n)
    enz = np.zeros(n)
    J4 = np.zeros((4, 4))
    J4[0, 1] = 1.0
    J4[1, 0] = -1.0
    J4[2, 3] = 1.0
    J4[3, 2] = -1.0
    for k in range(n):
        ev = np.linalg.eigvals(matrix[k, :4, :4] @ J4)
        imag = np.sort(np.abs(ev.imag))
        imag = imag[imag > 1e-30]
        if len(imag) >= 4:
            enx[k] = imag[0]
            eny[k] = imag[2]
        elif len(imag) >= 2:
            enx[k] = imag[0]
            eny[k] = imag[1]
        z2 = matrix[k, 4, 4]
        pz2 = matrix[k, 5, 5]
        zpz = matrix[k, 4, 5]
        enz[k] = float(np.sqrt(max(z2 * pz2 - zpz**2, 0.0)))

    name = Path(base).name
    logger.info("Read Sigma file '%s': %d rows", name, n)
    return SigmaData(z=z, e_kin_eV=e_kin, matrix=matrix,
                     enx=enx, eny=eny, enz=enz, filename=name)


def read_ref_file(rootname: str, run: str = "001") -> RefData:
    """Read a reference-particle trajectory file (9 columns)."""
    base = _base_from_rootname(rootname)
    path = Path(base + ".ref." + run)
    data = _load_columns(path, 9)
    name = Path(base).name
    return RefData(
        z=data[:, 0],
        t=data[:, 1] * NS_TO_S,
        pz=data[:, 2] * MEV_TO_EV,
        dedz=data[:, 3] * (MEV_TO_EV / 1.0),
        larmor=data[:, 4],
        xoff=data[:, 5] * MM_TO_M,
        yoff=data[:, 6] * MM_TO_M,
        px=data[:, 7],
        py=data[:, 8],
        filename=name,
    )


def read_log_file(rootname: str, run: str = "001") -> str:
    """Read an ASTRA log file as plain text."""
    base = _base_from_rootname(rootname)
    return Path(base + ".Log." + run).read_text()

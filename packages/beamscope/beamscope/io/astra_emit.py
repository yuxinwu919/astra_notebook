"""ASTRA Emit / Sigma / Ref / Log file readers.

All ASTRA evolution files are ASCII with space-separated scientific notation.
Read directly via np.loadtxt. Column layouts follow ASTRA Manual V3.2 Table 4.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger(__name__)


# ============================================================
# Data models
# ============================================================

@dataclass
class EmitData:
    """Single-plane emittance evolution data (7 columns).

    For Xemit/Yemit: avg=<x> or <y> [m], rms=sigma_x or sigma_y [m],
                      rmsprime=sigma_x' or sigma_y' [rad], emit=gamma*eps [m*rad]
    For Zemit: avg=E_kin [eV], rms=sigma_z [m], rmsprime=sigma_E [eV],
               emit=gamma*eps_z [eV*s]
    """
    z: np.ndarray          # [m] longitudinal position
    t: np.ndarray          # [ns] time of flight
    avg: np.ndarray        # <x>, <y>, or E_kin
    rms: np.ndarray        # sigma_x, sigma_y, or sigma_z
    rmsprime: np.ndarray   # sigma_x', sigma_y', or sigma_E
    emit: np.ndarray       # gamma * epsilon (normalized emittance)
    corr: np.ndarray       # correlation <x*x'> etc.
    label: str = ""        # 'x', 'y', or 'z'
    plane: str = ""        # 'horizontal', 'vertical', 'longitudinal'

    def to_dict(self) -> dict:
        return {
            "z": self.z, "t": self.t, "avg": self.avg,
            "rms": self.rms, "rmsprime": self.rmsprime,
            "emit": self.emit, "corr": self.corr,
            "label": self.label, "plane": self.plane,
        }


@dataclass
class EmitSet:
    """Full set of X/Y/Z emittance evolution + optional Cemit."""
    x: EmitData
    y: EmitData
    z: EmitData
    c: Optional[EmitData] = None
    filename: str = ""

    @property
    def has_cemit(self) -> bool:
        return self.c is not None


@dataclass
class SigmaData:
    """6x6 beam covariance matrix evolution (23 columns).

    Columns: z, gamma, then 21 upper-triangle elements of the symmetric
    6x6 covariance matrix: x2, xpx, xy, xpy, xz, xpz, px2, pxy, pxpy,
    pxz, pxpz, y2, ypy, yz, ypz, py2, pyz, pypz, z2, zpz, pz2
    """
    z: np.ndarray          # [m]
    gamma: np.ndarray      # relativistic gamma
    matrix: np.ndarray     # (N, 6, 6) full covariance
    enx: np.ndarray        # eigen-emittance x [m*rad]
    eny: np.ndarray        # eigen-emittance y [m*rad]
    enz: np.ndarray        # longitudinal emittance [m]
    filename: str = ""

    def to_dict(self) -> dict:
        return {
            "z": self.z, "gamma": self.gamma,
            "enx": self.enx, "eny": self.eny, "enz": self.enz,
            "filename": self.filename,
        }


@dataclass
class RefData:
    """Reference particle trajectory (9 columns).

    Columns follow ASTRA V3.2 Table 4:
    z[m], t[ns], pz[MeV/c], dE/dz[MeV/m], Larmor angle[rad],
    xoff[mm], yoff[mm], px[eV/c], py[eV/c]
    """
    z: np.ndarray          # [m]
    t: np.ndarray          # [ns]
    pz: np.ndarray         # [MeV/c]
    dE_dz: np.ndarray      # [MeV/m]
    larmor: np.ndarray     # [rad]
    xoff: np.ndarray       # [mm]
    yoff: np.ndarray       # [mm]
    px: np.ndarray         # [eV/c]
    py: np.ndarray         # [eV/c]
    filename: str = ""


@dataclass
class SimSet:
    """Complete set of simulation output files for one ASTRA run."""
    base_name: str
    emit: Optional[EmitSet] = None
    sigma: Optional[SigmaData] = None
    ref: Optional[RefData] = None
    log_text: Optional[str] = None
    phase_files: list[Path] = field(default_factory=list)
    dist_files: list[Path] = field(default_factory=list)


# ============================================================
# Readers
# ============================================================

def read_emit_files(rootname: str, run: str = "001") -> EmitSet:
    """Read Xemit, Yemit, Zemit files for a simulation run.

    Args:
        rootname: Base path (e.g., '/path/to/Example') or full path with stem.
        run: Run number string (e.g., '001').

    Returns:
        EmitSet with x, y, z emit data.
    """
    root = Path(rootname)
    base = str(root.parent / root.stem) if root.suffix else str(root)

    x_path = Path(f"{base}.Xemit.{run}")
    y_path = Path(f"{base}.Yemit.{run}")
    z_path = Path(f"{base}.Zemit.{run}")

    x_data = _read_emit_plane(x_path, "x", "horizontal")
    y_data = _read_emit_plane(y_path, "y", "vertical")
    z_data = _read_emit_plane(z_path, "z", "longitudinal")

    # Try Cemit
    c_data = None
    c_path = Path(f"{base}.Cemit.{run}")
    if c_path.exists():
        try:
            c_data = _read_cemit(c_path)
        except Exception:
            pass

    filename = root.stem if root.suffix else Path(base).name

    logger.info(
        f"Read emit files for '{filename}': "
        f"Xemit={len(x_data.z)} rows, Yemit={len(y_data.z)} rows, "
        f"Zemit={len(z_data.z)} rows"
    )

    return EmitSet(x=x_data, y=y_data, z=z_data, c=c_data, filename=filename)


def read_sigma_file(rootname: str, run: str = "001") -> SigmaData:
    """Read Sigma (6x6 beam covariance matrix) file.

    Args:
        rootname: Base path.
        run: Run number string.

    Returns:
        SigmaData with full covariance matrix and eigen-emittances.
    """
    root = Path(rootname)
    base = str(root.parent / root.stem) if root.suffix else str(root)
    path = Path(f"{base}.Sigma.{run}")

    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    z = data[:, 0]
    gamma = data[:, 1]
    cov_elements = data[:, 2:]  # (N, 21)

    n = len(z)
    matrix = np.zeros((n, 6, 6))

    # Upper-triangle indices for 6x6 (row-major)
    idx = 0
    for i in range(6):
        for j in range(i, 6):
            matrix[:, i, j] = cov_elements[:, idx]
            if i != j:
                matrix[:, j, i] = cov_elements[:, idx]
            idx += 1

    # Eigen-emittances from 4x4 transverse block (x, px, y, py)
    enx = np.zeros(n)
    eny = np.zeros(n)
    enz = np.zeros(n)

    for k in range(n):
        # 4x4 transverse block
        sigma_4x4 = matrix[k, :4, :4]
        # Symplectic unit matrix
        J4 = np.zeros((4, 4))
        J4[0, 1] = 1; J4[1, 0] = -1
        J4[2, 3] = 1; J4[3, 2] = -1

        eigenvalues = np.linalg.eigvals(sigma_4x4 @ J4)
        imag_parts = np.sort(np.abs(np.imag(eigenvalues)))
        # Take every other (pairs)
        if len(imag_parts) >= 4:
            enx[k] = imag_parts[0]
            eny[k] = imag_parts[2]
        elif len(imag_parts) >= 2:
            enx[k] = imag_parts[0]
            eny[k] = imag_parts[1]

        # Longitudinal emittance: sqrt(z2*pz2 - zpz^2)
        z2 = matrix[k, 4, 4]
        pz2 = matrix[k, 5, 5]
        zpz = matrix[k, 4, 5]
        enz[k] = np.sqrt(max(z2 * pz2 - zpz**2, 0))

    filename = root.stem if root.suffix else Path(base).name

    logger.info(
        f"Read sigma file '{filename}': {n} rows, "
        f"enx={enx[0]:.4e}, eny={eny[0]:.4e}, enz={enz[0]:.4e} (first row)"
    )

    return SigmaData(
        z=z, gamma=gamma, matrix=matrix,
        enx=enx, eny=eny, enz=enz,
        filename=filename,
    )


def read_ref_file(rootname: str, run: str = "001") -> RefData:
    """Read reference particle trajectory file.

    Columns per ASTRA V3.2 Table 4:
    z[m], t[ns], pz[MeV/c], dE/dz[MeV/m], Larmor angle[rad],
    xoff[mm], yoff[mm], px[eV/c], py[eV/c]
    """
    root = Path(rootname)
    base = str(root.parent / root.stem) if root.suffix else str(root)
    path = Path(f"{base}.ref.{run}")

    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    filename = root.stem if root.suffix else Path(base).name
    logger.info(f"Read ref file '{filename}': {len(data)} rows")

    return RefData(
        z=data[:, 0], t=data[:, 1], pz=data[:, 2],
        dE_dz=data[:, 3], larmor=data[:, 4],
        xoff=data[:, 5], yoff=data[:, 6],
        px=data[:, 7], py=data[:, 8],
        filename=filename,
    )


def read_log_file(rootname: str, run: str = "001") -> str:
    """Read ASTRA log file as plain text."""
    root = Path(rootname)
    base = str(root.parent / root.stem) if root.suffix else str(root)
    path = Path(f"{base}.Log.{run}")
    return path.read_text()


def discover_sim_files(sim_dir: Path) -> SimSet:
    """Auto-discover all ASTRA output files for a simulation run.

    Given a directory containing files like Example.ini, Example.Xemit.001, etc.,
    group them by base name and return a SimSet.

    Args:
        sim_dir: Directory containing ASTRA output files.

    Returns:
        SimSet with all discovered files. Use the first non-distribution base name found.
    """
    # Find all ASTRA output files
    all_files = sorted(sim_dir.glob("*"))

    # Group by base name (e.g., 'Example')
    bases: dict[str, dict] = {}
    for f in all_files:
        parts = f.name.split(".")
        if len(parts) >= 2:
            base = parts[0]
            ext = parts[1]
            if base not in bases:
                bases[base] = {"ini": [], "dist": [], "Xemit": None, "Yemit": None,
                               "Zemit": None, "Sigma": None, "ref": None, "Log": None}

            if ext == "ini":
                bases[base]["ini"].append(f)
            elif ext == "Xemit":
                bases[base]["Xemit"] = f
            elif ext == "Yemit":
                bases[base]["Yemit"] = f
            elif ext == "Zemit":
                bases[base]["Zemit"] = f
            elif ext == "Sigma":
                bases[base]["Sigma"] = f
            elif ext == "ref":
                bases[base]["ref"] = f
            elif ext == "Log":
                bases[base]["Log"] = f
            elif ext.isdigit():
                bases[base]["dist"].append(f)

    if not bases:
        raise FileNotFoundError(f"No ASTRA files found in {sim_dir}")

    # Take the first base with the most files
    best = max(bases.items(), key=lambda x: len([v for v in x[1].values() if v]))
    base_name, files = best

    sim_set = SimSet(base_name=base_name)
    sim_set.dist_files = files["dist"]

    # Read emit
    sim_path = sim_dir / base_name
    try:
        sim_set.emit = read_emit_files(str(sim_path), "001")
    except Exception as e:
        logger.warning(f"Could not read emit files for {base_name}: {e}")

    try:
        sim_set.sigma = read_sigma_file(str(sim_path), "001")
    except Exception as e:
        logger.warning(f"Could not read sigma file for {base_name}: {e}")

    try:
        sim_set.ref = read_ref_file(str(sim_path), "001")
    except Exception as e:
        logger.warning(f"Could not read ref file for {base_name}: {e}")

    try:
        sim_set.log_text = read_log_file(str(sim_path), "001")
    except Exception as e:
        logger.warning(f"Could not read log file for {base_name}: {e}")

    return sim_set


# ============================================================
# Internal helpers
# ============================================================

def _read_emit_plane(path: Path, label: str, plane: str) -> EmitData:
    """Read a single emit plane file (Xemit, Yemit, or Zemit).

    All three have the same 7-column format:
    z[m], t[ns], avg, rms, rmsprime, emit, corr
    """
    if not path.exists():
        raise FileNotFoundError(f"Emit file not found: {path}")

    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    return EmitData(
        z=data[:, 0], t=data[:, 1],
        avg=data[:, 2], rms=data[:, 3],
        rmsprime=data[:, 4], emit=data[:, 5], corr=data[:, 6],
        label=label, plane=plane,
    )


def _read_cemit(path: Path) -> EmitData:
    """Read Cemit (core emittance) file. 13 columns: z + 12 emittance values."""
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)

    # This is a simplified representation; full Cemit support would
    # store all 12 core emittance columns separately.
    return EmitData(
        z=data[:, 0], t=np.zeros(len(data)),
        avg=data[:, 1],  # emit100x as placeholder
        rms=data[:, 5],  # emit100y as placeholder
        rmsprime=data[:, 9],  # emit100z as placeholder
        emit=data[:, 1],  # emit100x
        corr=data[:, 5],  # emit100y
        label="core", plane="all",
    )

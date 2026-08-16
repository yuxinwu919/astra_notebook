"""ASTRA evolution-file readers: Xemit/Yemit/Zemit/Cemit/Sigma/ref/Log.

All formats per ASTRA Manual V3.2, Table 4. Units are converted to SI
on read; the original display units are documented per column.

    Xemit / Yemit   (7 cols): z[m] t[ns] u_avr[mm] u_rms[mm]
                              u'_rms[mrad] eps_n[1e-6 m.rad] <uu'>_avr/sigma_u[mrad]
    Zemit           (7 cols): z[m] t[ns] E_kin[MeV] z_rms[mm]
                              dE_rms[keV] eps_zn[keV.mm] <z E'>_avr[keV]
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
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

import numpy as np

from ..constants import M_E_C2_EV

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
    emit=eps_n [m.rad], corr=<u u'>/sigma_u [rad] (手册第 7 列语义).
    For Zemit: avg=E_kin [eV], rms=sigma_z [m], rmsprime=sigma_E [eV],
    emit=eps_zn [eV.m], corr=<z E'> [eV].
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
    manual 4.13. 文件把动量列与能量列都归一化到 mc
    (无量纲); 方差元 (能量平方) 随之是 mc^2 因子 — 历史"3.83 因子"
    即 1/mc^2 = 3.83 (1/mc[MeV]^2);
    读取时统一转换为 SI: 位置 [m], 动量 [eV/c], 能量 [eV]。
    """

    z: np.ndarray          # [m]
    e_kin_eV: np.ndarray   # [eV] (file stores MeV)
    matrix: np.ndarray     # (N, 6, 6) covariance, SI
    enx: np.ndarray        # 归一化 eigen-emittance x [m.rad] (可与 Xemit 对照)
    eny: np.ndarray        # 归一化 eigen-emittance y [m.rad]
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
        corr = data[:, 6] * KEV_TO_EV  # <z E'>_avr [keV -> eV] (manual Table 4)
    else:
        avg = data[:, 2] * MM_TO_M          # <u> [m]
        rms = data[:, 3] * MM_TO_M          # sigma_u [m]
        rmsprime = data[:, 4] * MRAD_TO_RAD  # sigma_u' [rad]
        emit = data[:, 5] * 1e-6            # eps_n [m.rad] (no pi factor!)
        corr = data[:, 6] * MRAD_TO_RAD     # cov(u,u')/sigma_u [rad]

    return EmitData(z=z, t=t, avg=avg, rms=rms, rmsprime=rmsprime,
                    emit=emit, corr=corr, label=label, plane=plane)


def read_cemit_file(path) -> dict:
    """Read a Cemit (core emittance) file: 13 columns (批 3 统一).

    与 parse_output_file 同一表驱动实现, 返回标准化 dict:
      mean_z, norm_emit_x / core_emit_95percent_x / ... / norm_emit_z ...
    横向 m.rad (文件值 x1e-6), 纵向 eV.m (文件值 x1, keV.mm 数值一致)。
    此前这里用 EmitData + 动态 _cemit 属性塞 8 列 (数据模型 hack),
    与表驱动读取器并存且单位表重复 — 已收敛为单一实现。
    """
    return parse_output_file(path)


def read_sigma_file(rootname: str, run: str = "001") -> SigmaData:
    """Read a Sigma (6x6 covariance) file: 23 columns.

    z[m], E_kin[MeV], then the 21 upper-triangle elements sig(i,j),
    i=1..6, j=i..6 (row-major) of the covariance matrix in canonical
    coordinates (x, p~x, y, p~y, z, E_kin) - manual 4.13.

    验证 (examples/Manual_Example, 全部行):
      * sig(1,1) [m^2] = sigma_x^2, sig(5,5) = sigma_z^2 (与
        Xemit/Zemit 精确一致)
      * sig(2,2)/(mc)^2 = (sigma_x' * p_ref)^2, sig(6,6)/(mc)^2 =
        sigma_E^2 - 文件把动量列与能量列都归一化到 mc (方差元 mc^2);
        历史上"3.83 因子"即 1/mc^2 = 3.83 (见 physics_notes/06)
      * 由 Sigma 导出的归一化 eigen-emittance 与 Xemit 的 eps_n
        逐行对照 < 2% (test/test_cross_validation.py)
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

    # 文件单位 -> SI: 列 2/4 (p~x/p~y) 与列 6 (E_kin) 均"除以 mc"。
    # 协方差按元素对缩放: matrix[i,j] *= s_i * s_j, 因此对角元
    # sig66 乘 mc^2, 混合元 (z,E) 乘 mc。
    mc = M_E_C2_EV
    scale = np.array([1.0, mc, 1.0, mc, 1.0, mc])
    matrix = matrix * scale[None, :, None] * scale[None, None, :]

    # Eigen-emittances: imaginary parts of eigenvalues of Sigma @ J
    # (4x4 transverse block, 单位 m.eV/c)。除以 mc 得到归一化
    # 发射度 [m.rad] (beta*gamma = p/mc 已含在换算中)。
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
            enx[k] = imag[0] / mc
            eny[k] = imag[2] / mc
        elif len(imag) >= 2:
            enx[k] = imag[0] / mc
            eny[k] = imag[1] / mc
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
    return Path(base + ".Log." + run).read_text(encoding="utf-8", errors="replace")


# ============================================================
# Generic tabular output parser
# ============================================================
# Column name/factor/unit tables vendored and adapted from lume-astra
# (C. Mayes et al., LBNL, Apache-2.0):
#   https://github.com/ChristopherMayes/lume-astra  (astra/parsers.py)
# Original units = what ASTRA writes (manual Table 4); factors convert
# to SI. Emittance columns: 'mm-mrad'/'mm-keV' with factor 1e-6 / 1,
# numerically equal to the pi-unit values printed by ASTRA.

OUTPUT_TABLES = {
    "Xemit": (
        ["mean_z", "mean_t", "mean_x", "sigma_x", "sigma_xp",
         "norm_emit_x", "cov_x__xp/sigma_x"],
        ["m", "ns", "mm", "mm", "mrad", "mm-mrad", "mrad"],
        [1, 1e-9, 1e-3, 1e-3, 1e-3, 1e-6, 1e-3],
        ["m", "s", "m", "m", "1", "m", "rad"],
    ),
    "Yemit": (
        ["mean_z", "mean_t", "mean_y", "sigma_y", "sigma_yp",
         "norm_emit_y", "cov_y__yp/sigma_y"],
        ["m", "ns", "mm", "mm", "mrad", "mm-mrad", "mrad"],
        [1, 1e-9, 1e-3, 1e-3, 1e-3, 1e-6, 1e-3],
        ["m", "s", "m", "m", "1", "m", "rad"],
    ),
    "Zemit": (
        ["mean_z", "mean_t", "mean_kinetic_energy", "sigma_z",
         "sigma_energy", "norm_emit_z", "cov_z__energy/sigma_z"],
        ["m", "ns", "MeV", "mm", "keV", "mm-keV", "keV"],
        [1, 1e-9, 1e6, 1e-3, 1e3, 1, 1e3],
        ["m", "s", "eV", "m", "eV", "m*eV", "eV"],
    ),
    "Cemit": (
        ["mean_z", "norm_emit_x", "core_emit_95percent_x",
         "core_emit_90percent_x", "core_emit_80percent_x",
         "norm_emit_y", "core_emit_95percent_y",
         "core_emit_90percent_y", "core_emit_80percent_y",
         "norm_emit_z", "core_emit_95percent_z",
         "core_emit_90percent_z", "core_emit_80percent_z"],
        ["m"] + 8 * ["mm-mrad"] + 4 * ["mm-keV"],
        # NOTE: deviates from lume-astra, which uses 1e-6 for ALL 12
        # columns. The z-plane is 'keV mm' = 1 eV.m, so the factor is 1
        # (lume-astra's copy-paste bug would make eps_zn 1e6 times too
        # small). Verified against ASTRA's own Cemit.001 values.
        [1] + 8 * [1e-6] + 4 * [1],
        ["m"] + 8 * ["m"] + 4 * ["m*eV"],
    ),
    "LandF": (
        ["landf_z", "landf_n_particles", "landf_total_charge",
         "landf_n_lost", "landf_energy_deposited", "landf_energy_exchange"],
        ["m", "1", "nC", "1", "J", "J"],
        [1, 1, 1e-9, 1, 1, 1],
        ["m", "1", "C", "1", "J", "J"],
    ),
}

# suffix (lower) -> table key
_EXT_TO_TABLE = {
    ".xemit": "Xemit", ".yemit": "Yemit", ".zemit": "Zemit",
    ".cemit": "Cemit", ".landf": "LandF",
}


def output_file_type(path) -> str:
    """ASTRA output type from filename: 'Example.Xemit.001' -> 'Xemit'."""
    parts = str(Path(path)).rsplit(".", 2)
    if len(parts) == 3 and parts[2].isdigit():
        return parts[1]
    return parts[-1]


def parse_output_file(path, standardize_labels: bool = True) -> dict:
    """Parse any tabular ASTRA output file into SI unit arrays.

    Covers Xemit/Yemit/Zemit/Cemit/LandF (Table 4). Returns a dict of
    1D arrays keyed by column name. With standardize_labels, the
    covariance columns (cov_x__xp/sigma_x etc.) are multiplied back to
    plain covariances: cov_x__xp, cov_y__yp, cov_z__energy.
    """
    path = Path(path)
    data = np.loadtxt(path, ndmin=2)
    if data.size == 0:
        raise ValueError("empty output file: " + str(path))

    ftype = output_file_type(path)
    key = _EXT_TO_TABLE.get("." + ftype.lower())
    if key is None:
        raise ValueError("unsupported output type: " + ftype)

    names, _, factors, _ = OUTPUT_TABLES[key]
    if data.shape[1] < len(names):
        # 批 3: 截断文件显式报错 (此前静默返回部分列)
        raise ValueError(
            "%s: 需要 >= %d 列, 实际 %d"
            % (path, len(names), data.shape[1]))
    d = {}
    for i, name in enumerate(names):
        d[name] = data[:, i] * factors[i]

    if standardize_labels:
        if ftype == "Xemit" and "cov_x__xp/sigma_x" in d:
            d["cov_x__xp"] = d.pop("cov_x__xp/sigma_x") * d["sigma_x"]
        if ftype == "Yemit" and "cov_y__yp/sigma_y" in d:
            d["cov_y__yp"] = d.pop("cov_y__yp/sigma_y") * d["sigma_y"]
        if ftype == "Zemit" and "cov_z__energy/sigma_z" in d:
            d["cov_z__energy"] = d.pop("cov_z__energy/sigma_z") * d["sigma_z"]
    return d

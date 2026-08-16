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

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..constants import C_LIGHT


@dataclass
class CavityField:
    """On-axis cavity field table.

    z [m], ez0 原始任意单位 (手册 6.9: 场表第 2 列为任意单位,
    峰值按 deck 的 MaxE 缩放; C_noscale=T 时数值才直接是 MV/m)。
    """

    z: np.ndarray
    ez0: np.ndarray
    source: str = ""

    @property
    def peak_field_MVpm(self) -> float:
        return float(np.max(np.abs(self.ez0)))   # 任意单位峰值

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
    blocks: list = field(default_factory=list)   # 多块读取时的全部块 (批 5: 声明字段)


def read_cavity_field(path) -> CavityField:
    """Read a two-column cavity field table (z [m], Ez [MV/m])."""
    path = Path(path)
    data = np.loadtxt(path)
    if data.ndim == 1:
        data = data.reshape(1, -1)
    if data.shape[1] < 2:
        raise ValueError("cavity field file must have >= 2 columns: " + str(path))
    return CavityField(z=data[:, 0].astype(float),
                       ez0=data[:, 1].astype(float),          # 原始任意单位 (手册 6.9)
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
    """Read an ASTRA wake table (manual section 6.8).

    两种格式 (批 3 起都支持):
      多块: 首行 (nblocks, 0), 随后每块一行头 (N, 0) 加 N 行 (s, W);
      单块 (手册 6.8 原文): 首行 (N, 0) 即块头, 后接 N 行 (s, W)。
    返回第一块, 所有块挂在其 .blocks 属性上。
    """
    path = Path(path)
    lines = [ln.split() for ln in path.read_text(encoding="utf-8").splitlines()
             if ln.strip() and not ln.strip().startswith(("#", "!"))]
    if not lines:
        raise ValueError("empty wake file: " + str(path))
    if len(lines[0]) < 2 or float(lines[0][1]) != 0:
        raise ValueError("wake header must be '<N> 0': " + str(path))
    nblocks = int(float(lines[0][0]))
    if nblocks <= 0:
        raise ValueError("wake file declares no blocks: " + str(path))

    # 单块探测: 第二行不是 "(N, 0)" 块头 -> 首行即单块头 "N 0"
    single = (
        len(lines) < 2
        or len(lines[1]) != 2
        or float(lines[1][1]) != 0
    )
    blocks: list = []
    if single:
        n = nblocks
        rows_raw = lines[1:1 + n]
        if len(rows_raw) < n:
            raise ValueError("wake block truncated: " + str(path))
        rows = np.asarray([[float(v) for v in ln] for ln in rows_raw],
                          dtype=float)
        if rows.ndim != 2 or rows.shape[1] < 2:
            raise ValueError("wake block truncated: " + str(path))
        blocks.append(WakePotential(s=rows[:, 0], w=rows[:, 1],
                                     source=str(path)))
    else:
        i = 1
        for _ in range(nblocks):
            if i >= len(lines):
                raise ValueError(
                    "wake file truncated (block header): " + str(path))
            n = int(float(lines[i][0]))
            i += 1
            rows = np.asarray([[float(v) for v in ln] for ln in lines[i:i + n]],
                              dtype=float)
            if rows.ndim != 2 or rows.shape[0] < n or rows.shape[1] < 2:
                raise ValueError("wake block truncated: %s" % path)
            blocks.append(WakePotential(s=rows[:, 0], w=rows[:, 1],
                                         source=str(path)))
            i += n
    out = blocks[0]
    out.blocks = blocks
    return out


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
    with open(path, encoding="utf-8") as f:
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
    """修复 MATLAB 写的 laser.dat 3D 图头 (就地转换)。

    DESY Plasma_Example_2 的 laser.dat 头三行为 (n, min, spacing),
    但 n 写成浮点形式 (8.1e+01)。实测 (macOS Apple Silicon ASTRA
    构建):
      * 浮点计数 8.1e+01            -> "Error while reading file"
      * 逐值网格行 (n 后跟 n 个值)   -> 读取 3D 图时 SIGSEGV
      * 整数计数 + (n,min,spacing)  -> 正常 (本函数的输出)

    因此只把计数转成整数, 其余原样保留。
    """
    path = Path(path)
    with open(path, encoding="utf-8") as f:
        lines = f.readlines()
    if len(lines) < 4:
        raise ValueError("laser map too short: " + str(path))
    out = []
    for i in range(3):
        vals = lines[i].split()
        n = int(float(vals[0]))
        out.append("%d %s %s\n" % (n, vals[1], vals[2]))
    with open(path, "w", encoding="utf-8") as f:
        f.writelines(out + lines[3:])
    return out

def read_3d_field_map(path):
    """读取 ASTRA 3D 场图 (如 3D_test.ex / 3D_Dipole.bx / laser.dat).

    网格头支持两种形式 (每行一个网格, 共 3 行):
      * 逐值: n 后跟 n 个网格值 (手册格式)
      * 紧凑: n, min, spacing (MATLAB/DESY 常见, 就地展开)
    随后为数据值 (可跨行), 顺序 x 最快、y 次之、z 最慢 (Fortran 序)。

    Returns:
        (x, y, z, F) — F 为 (nx, ny, nz) 数组, 按 F[ix, iy, iz] 索引,
        SI 单位按文件名约定 (ex/ey/ez: V/m; bx/by/bz: T; 数值原样返回)。
    """
    path = Path(path)
    # 批 5: 单次读入, 全局 token 流 (此前紧凑分支读一遍、逐值分支
    # 又 read_text() 读第二遍并逐 token 转 float, 65MB laser.dat 翻倍)
    vals = [float(v) for v in path.read_text(encoding="utf-8").split()]
    if len(vals) < 4:
        raise ValueError("3D map too short: " + str(path))

    # 紧凑头 (n, min, spacing) 候选: 前 9 个 token 恰为
    # n1,min1,dx1, n2,min2,dx2, n3,min3,dx3 且各自 n>=1,
    # 且剩余数据长度与网格积精确相等 — 才按紧凑头解析;
    # 否则按逐值头解析 (3 个 token 可能只是恰好每行 3 个网格值
    # 的巧合, 如 3D_Dipole.bx)
    head9 = vals[:9]
    compact = False
    if all(float(h) >= 1 for h in (head9[0], head9[3], head9[6])):
        ns = [int(head9[0]), int(head9[3]), int(head9[6])]
        need = ns[0] * ns[1] * ns[2]
        if len(vals) - 9 == need:
            compact = True
            grids = [np.array([head9[1] + j * head9[2] for j in range(ns[0])],
                              dtype=float),
                     np.array([head9[4] + j * head9[5] for j in range(ns[1])],
                              dtype=float),
                     np.array([head9[7] + j * head9[8] for j in range(ns[2])],
                              dtype=float)]
            data = np.array(vals[9:], dtype=float)
    if not compact:
        # 逐值头: 全局 token 流, 支持换行 (自由格式)
        idx = 0
        grids = []
        for _ in range(3):
            n = int(vals[idx])
            if n < 1:
                raise ValueError("3D map grid size invalid: " + str(path))
            idx += 1
            grids.append(np.array(vals[idx:idx + n], dtype=float))
            idx += n
        data = np.array(vals[idx:], dtype=float)
    x, y, z = grids
    nx, ny, nz = len(x), len(y), len(z)
    if len(data) < nx * ny * nz:
        raise ValueError("3D map data truncated: " + str(path))
    f = data[:nx * ny * nz].reshape(nx, ny, nz, order="F")
    return x, y, z, f

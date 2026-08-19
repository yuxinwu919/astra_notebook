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

Off-axis field expansion (manual V3.2, Appendix I; 3rd order like
ASTRA's C_higher_order/S_higher_order=TRUE default):
    TM (cavity, amplitude; sin(wt) factors omitted):
        Ez(r,z)  = Ez0 - (r^2/4)(Ez0'' + w^2/c^2 Ez0)
        Er(r,z)  = -(r/2) Ez0' + (r^3/16)(Ez0''' + w^2/c^2 Ez0')
        Bphi(r,z)= [(r/2) Ez0 - (r^3/16)(Ez0'' + w^2/c^2 Ez0)] w/c^2
    solenoid (3rd order):
        Bz(r,z)  = Bz0 - (r^2/4) Bz0'' + (r^4/64) Bz0''''
        Br(r,z)  = -(r/2) Bz0' + (r^3/16) Bz0'''
    (2026-08 audit: the old TM branch used a static expansion without
    the w^2/c^2 terms and with a spurious r^4 term - 19.6% Ez error at
    r=2 cm on the real 3_cell_L-Band.dat.)
"""

from __future__ import annotations

import warnings
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
    def peak_field_arb(self) -> float:
        """任意单位峰值 (批 5: 不再伪装 MV/m, 见手册 6.9)。"""
        return float(np.max(np.abs(self.ez0)))

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
        ez0p = np.interp(z, self.z, dez)
        ez0pp = np.interp(z, self.z, d2ez)
        ez0ppp = np.interp(z, self.z, d3ez)
        # 手册附录 I (TM 驻波, 3 阶; 2026-08 审计 P1: 旧式为静态展开,
        # 缺 w²/c² 修正且多 r⁴ 项)。静场 (omega=0) 时 w²/c² 项自然为零。
        w2c2 = (omega / C_LIGHT) ** 2 if omega else 0.0
        ez = ez0 - (r**2 / 4.0) * (ez0pp + w2c2 * ez0)
        er = -(r / 2.0) * ez0p + (r**3 / 16.0) * (ez0ppp + w2c2 * ez0p)
        if omega:
            bphi = ((r / 2.0) * ez0
                    - (r**3 / 16.0) * (ez0pp + w2c2 * ez0)) * (omega / C_LIGHT**2)
        else:
            bphi = np.zeros_like(r)
        return ez, er, bphi

    def expansion_radius(self, omega: float = 0.0, smooth_window=None):
        """场展开半径 R3rd (手册 8 章, TM) vs z.

        R3rd_Er = sqrt(|0.08 Ez'| / |Ez''' + w²/c² Ez'|)
        R3rd_Bφ = sqrt(|0.08 Ez0| / |Ez'' + w²/c² Ez0|)

        smooth_window: 可选 Savitzky-Golay 平滑窗口 (奇数, 手册 8 章
        C_smooth 的演示), 抑制场表数值噪声导致的 R3rd 尖峰; 默认
        None = 不平滑 (纯数值)。

        返回 (r3rd_er, r3rd_bphi) [m]; 数值噪声处可能 inf/NaN。
        """
        ez0 = self._smooth(self.ez0, smooth_window)
        dez = np.gradient(ez0, self.z)
        d2ez = np.gradient(dez, self.z)
        d3ez = np.gradient(d2ez, self.z)
        w2c2 = (omega / C_LIGHT) ** 2
        with np.errstate(divide="ignore", invalid="ignore"):
            r_er = np.sqrt(np.abs(0.08 * dez) / np.abs(d3ez + w2c2 * dez))
            r_bp = np.sqrt(np.abs(0.08 * ez0)
                           / np.abs(d2ez + w2c2 * ez0))
        return r_er, r_bp

    @staticmethod
    def _smooth(v, window):
        """可选 Savitzky-Golay 平滑 (奇数窗口, polyorder=3)."""
        if not window or len(v) < window or window < 3:
            return v
        from scipy.signal import savgol_filter
        return savgol_filter(v, int(window) | 1, 3)


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
        if peak == 0.0:
            raise ValueError("cannot scale a zero-peak solenoid field")
        out = SolenoidField(z=self.z.copy(), bz0=self.bz0 * (max_b_T / peak),
                            source=self.source)
        return out

    def field_at(self, r: np.ndarray, z: np.ndarray):
        """Off-axis (Br, Bz) via the axis expansion, 3rd order
        (ASTRA S_higher_order=TRUE 默认; 手册附录 I).

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
        d3bz = np.gradient(d2bz, self.z)
        d4bz = np.gradient(d3bz, self.z)
        bz0p = np.interp(z, self.z, dbz)
        bz0pp = np.interp(z, self.z, d2bz)
        bz0ppp = np.interp(z, self.z, d3bz)
        bz0pppp = np.interp(z, self.z, d4bz)
        # 3 阶 (2026-08 审计 P3: 旧实现仅 1 阶, 与 ASTRA 默认不符)
        br = -(r / 2.0) * bz0p + (r**3 / 16.0) * bz0ppp
        bz = bz0 - (r**2 / 4.0) * bz0pp + (r**4 / 64.0) * bz0pppp
        return br, bz

    def expansion_radius(self, smooth_window=None):
        """R3rd = sqrt(|0.08 Bz'| / |Bz'''|) (手册 8 章, 静磁) [m].

        smooth_window: 可选 Savitzky-Golay 平滑 (抑制 Bz''' 噪声尖峰),
        默认 None = 不平滑。
        """
        bz0 = self._smooth(self.bz0, smooth_window)
        dbz = np.gradient(bz0, self.z)
        d3bz = np.gradient(np.gradient(dbz, self.z), self.z)
        with np.errstate(divide="ignore", invalid="ignore"):
            r = np.sqrt(np.abs(0.08 * dbz) / np.abs(d3bz))
        return r

    @staticmethod
    def _smooth(v, window):
        """可选 Savitzky-Golay 平滑 (奇数窗口, polyorder=3)."""
        if not window or len(v) < window or window < 3:
            return v
        from scipy.signal import savgol_filter
        return savgol_filter(v, int(window) | 1, 3)


@dataclass
class WakePotential:
    """Wake potential table: s [m], W [V/C] (pseudo-Green function)."""

    s: np.ndarray
    w: np.ndarray
    source: str = ""
    blocks: list = field(default_factory=list)   # 多块读取时的全部块 (批 5: 声明字段)


def read_cavity_field(path) -> CavityField:
    """Read a two-column cavity field table (z [m], Ez [arb. units]).

    手册 6.9: 场表第 2 列为任意单位, 峰值按 deck 的 MaxE 缩放;
    C_noscale=T 时数值才直接是 MV/m。z 必须严格单调递增。
    """
    path = Path(path)
    data = np.loadtxt(path, ndmin=2)
    if data.shape[1] < 2:
        raise ValueError("cavity field file must have >= 2 columns: " + str(path))
    zcol = data[:, 0].astype(float)
    if not np.all(np.diff(zcol) > 0):
        raise ValueError(
            "cavity field z must be strictly increasing: " + str(path)
            + " (2026-08 audit: np.interp/gradient silently return garbage "
            "on non-monotonic tables)")
    return CavityField(z=zcol,
                       ez0=data[:, 1].astype(float),          # 原始任意单位 (手册 6.9)
                       source=str(path))


def read_solenoid_field(path) -> SolenoidField:
    """Read a two-column solenoid field table (z [m], Bz [arb. units])."""
    path = Path(path)
    data = np.loadtxt(path, ndmin=2)
    if data.shape[1] < 2:
        raise ValueError("solenoid field file must have >= 2 columns: " + str(path))
    return SolenoidField(z=data[:, 0].astype(float),
                         bz0=data[:, 1].astype(float),
                         source=str(path))


@dataclass
class TEField:
    """TE 模腔场 (手册 6.9: 文件名以 'TE_' 开头, 表存轴上纵向磁场 Bz).

    z [m], bz0 任意单位; 按 deck 的 MaxE 缩放 (TE 模下 MaxE 指轴上
    纵向磁场分量)。离轴场按手册 8 章 TE 展开 (含 w²/c² 项)。
    """

    z: np.ndarray
    bz0: np.ndarray
    source: str = ""

    @property
    def peak_field_arb(self) -> float:
        return float(np.max(np.abs(self.bz0)))

    def scaled(self, max_field) -> "TEField":
        """缩放使峰值等于 max_field (TE 模的 MaxE 指轴上 Bz)."""
        peak = float(np.max(np.abs(self.bz0)))
        if peak == 0.0:
            raise ValueError("cannot scale a zero-peak TE field")
        return TEField(z=self.z.copy(), bz0=self.bz0 * (max_field / peak),
                       source=self.source)

    def field_at(self, r: np.ndarray, z: np.ndarray, omega: float = 0.0):
        """离轴 TE 展开 (手册 8 章): (bz, br, ephi) [T], [T], [V/m].

        Bz(r) = Bz0 - (r²/4)(Bz0'' + w²/c² Bz0)
        Br(r) = -(r/2) Bz0' + (r³/16)(Bz0''' + w²/c² Bz0')
        Eφ(r) = [(r/2) Bz0 - (r³/16)(Bz0'' + w²/c² Bz0)] ω
        """
        r = np.asarray(r, dtype=float)
        z = np.asarray(z, dtype=float)
        bz0 = np.interp(z, self.z, self.bz0)
        dbz = np.gradient(self.bz0, self.z)
        d2bz = np.gradient(dbz, self.z)
        d3bz = np.gradient(d2bz, self.z)
        bz0p = np.interp(z, self.z, dbz)
        bz0pp = np.interp(z, self.z, d2bz)
        bz0ppp = np.interp(z, self.z, d3bz)
        w2c2 = (omega / C_LIGHT) ** 2
        bz = bz0 - (r**2 / 4.0) * (bz0pp + w2c2 * bz0)
        br = -(r / 2.0) * bz0p + (r**3 / 16.0) * (bz0ppp + w2c2 * bz0p)
        ephi = ((r / 2.0) * bz0 - (r**3 / 16.0) * (bz0pp + w2c2 * bz0)) * omega
        return bz, br, ephi

    def expansion_radius(self, omega: float = 0.0, smooth_window=None):
        """R3rd for TE (手册 8 章): (r3rd_ephi, r3rd_br) [m].

        smooth_window: 可选 Savitzky-Golay 平滑, 默认 None = 不平滑。
        """
        bz0 = self._smooth(self.bz0, smooth_window)
        dbz = np.gradient(bz0, self.z)
        d2bz = np.gradient(dbz, self.z)
        w2c2 = (omega / C_LIGHT) ** 2
        with np.errstate(divide="ignore", invalid="ignore"):
            r_ep = np.sqrt(np.abs(0.08 * bz0)
                           / np.abs(d2bz + w2c2 * bz0))
            r_br = np.sqrt(np.abs(0.08 * dbz) / np.abs(d2bz + w2c2 * dbz))
        return r_ep, r_br

    @staticmethod
    def _smooth(v, window):
        """可选 Savitzky-Golay 平滑 (奇数窗口, polyorder=3)."""
        if not window or len(v) < window or window < 3:
            return v
        from scipy.signal import savgol_filter
        return savgol_filter(v, int(window) | 1, 3)


def read_te_field(path) -> TEField:
    """Read a two-column TE-mode field table (z [m], Bz [arb.]), 手册 6.9."""
    path = Path(path)
    data = np.loadtxt(path, ndmin=2)
    if data.shape[1] < 2:
        raise ValueError("TE field file must have >= 2 columns: " + str(path))
    return TEField(z=data[:, 0].astype(float),
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

    def _parse_multi():
        blocks = []
        i = 1
        for _ in range(nblocks):
            if i >= len(lines):
                raise ValueError(
                    "wake file truncated (block header): " + str(path))
            n = int(float(lines[i][0]))
            if len(lines[i]) < 2 or float(lines[i][1]) != 0 or n < 1:
                raise ValueError("wake block header invalid: " + str(path))
            i += 1
            rows = np.asarray([[float(v) for v in ln] for ln in lines[i:i + n]],
                              dtype=float)
            if rows.ndim != 2 or rows.shape[0] < n or rows.shape[1] < 2:
                raise ValueError("wake block truncated: %s" % path)
            blocks.append(WakePotential(s=rows[:, 0], w=rows[:, 1],
                                         source=str(path)))
            i += n
        if i != len(lines):
            raise ValueError("wake file has trailing data: %s" % path)
        return blocks

    def _parse_single():
        n = nblocks
        rows_raw = lines[1:1 + n]
        if len(rows_raw) < n or len(lines) != 1 + n:
            raise ValueError("wake block truncated: " + str(path))
        rows = np.asarray([[float(v) for v in ln] for ln in rows_raw],
                          dtype=float)
        if rows.ndim != 2 or rows.shape[1] < 2:
            raise ValueError("wake block truncated: " + str(path))
        return [WakePotential(s=rows[:, 0], w=rows[:, 1],
                              source=str(path))]

    # 批 5 复核修正: 先尝试多块解析 (含 nblocks=1), 失败再按单块;
    # 两者都能解析时以多块为准 (单块文件首数据行 W=0 会被旧探测误判)
    try:
        blocks = _parse_multi()
    except ValueError:
        blocks = _parse_single()
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
    """把含 m 个单元的 TWS 周期段 (z1..z2) 周期展开到 C_numb=n_cell 个单元。

    手册 6.9 (TWS): 场表头 4 值 (z1 z2 n m), z1..z2 为含 m 个单元的
    周期段; 按单元平移拼接, 布局

        | Entrance | Cells × (n_cell / m) | Exit |

    入口段 (zmin..z1) 与出口段 (z2..zmax) 的原表采样原样保留 (不动),
    仅周期段重复; 周期段在原表 [z1, z2) 内按原表密度重采样, 场值用
    np.interp 从原表插值 (原表必须覆盖一个完整周期段)。

    仅轴上场演示: 一阶横向展开不在本函数范围 (手册 6.9 要求
    C_numb·n/m 为整数; 与 2026-08 审计结论一致 — 本函数是示例绘图
    用的表级展开工具, 不是物理场展开的替代; n<0 返波不支持)。

    Args:
        z0: 原表 z [m], 严格递增, 覆盖完整周期段。
        f0: 原表轴上场值 (与 z0 等长)。
        z1, z2: 周期段起止 [m] (含 m_cells_in_body 个单元)。
        m_cells_in_body: 周期段内单元数 (头文件第 4 值 m)。
        n_cell: 总单元数 C_numb。

    Returns:
        (zfull, ffull) — 展开后的 (z, f) 表, 布局 |Entrance|
        Cells×(n_cell/m) |Exit|。

    Raises:
        ValueError: n_cell % m_cells_in_body != 0 (手册整数条件,
            不再静默截断), 或输入不合法 (z0 非递增 / 周期段越界 /
            表未覆盖周期段 / z2 <= z1)。
    """
    z0 = np.asarray(z0, dtype=float)
    f0 = np.asarray(f0, dtype=float)
    if z0.ndim != 1 or f0.ndim != 1 or z0.shape != f0.shape or len(z0) < 2:
        raise ValueError("z0/f0 必须为等长一维数组 (>= 2 点)")
    if not np.all(np.diff(z0) > 0):
        raise ValueError(
            "z0 必须严格递增 (2026-08 审计: 非单调表经 np.interp 会静默出垃圾)")
    if m_cells_in_body < 1:
        raise ValueError("m_cells_in_body 必须 >= 1, 得到 %s"
                         % (m_cells_in_body,))
    if n_cell < 1:
        raise ValueError("n_cell 必须 >= 1, 得到 %s" % (n_cell,))
    if n_cell % m_cells_in_body != 0:
        raise ValueError(
            "n_cell=%s 必须能被 m_cells_in_body=%s 整除 (手册 6.9: "
            "C_numb·n/m 需为整数); 不再静默截断"
            % (n_cell, m_cells_in_body))
    zmin, zmax = float(z0[0]), float(z0[-1])
    if z2 <= z1:
        raise ValueError("周期段长度 z2-z1 必须 > 0, 得到 %g" % (z2 - z1))
    if not (zmin <= z1 and z2 <= zmax):
        raise ValueError(
            "周期段 z1..z2 = [%g, %g] 必须在场表范围 [%g, %g] 内"
            % (z1, z2, zmin, zmax))
    n_body_samples = int(np.count_nonzero((z0 >= z1) & (z0 < z2)))
    if n_body_samples < 2:
        raise ValueError("场表必须覆盖一个完整周期段 (z1..z2 内至少 2 个采样点)")

    n_repeat = int(n_cell) // int(m_cells_in_body)
    l_period = z2 - z1

    # 入口/出口段原样保留 (手册 6.9: 仅周期段重复)
    mask_entrance = z0 < z1
    mask_exit = z0 >= z2
    z_entrance, f_entrance = z0[mask_entrance], f0[mask_entrance]
    z_exit, f_exit = z0[mask_exit], f0[mask_exit]

    # 周期段: 按原表采样密度在 [z1, z2) 建均匀网格, 场值 np.interp
    # 从原表取 (表格覆盖完整周期段, 展开才物理有效)
    z_body = np.linspace(z1, z2, n_body_samples, endpoint=False)
    f_body = np.interp(z_body, z0, f0)

    ztot = [z_entrance]
    ftot = [f_entrance]
    for i in range(n_repeat):
        ztot.append(z_body + i * l_period)
        ftot.append(f_body)
    ztot.append(z_exit + (n_repeat - 1) * l_period)
    ftot.append(f_exit)
    return np.concatenate(ztot), np.concatenate(ftot)


def fix_laser_map_header(path):
    """修复 MATLAB 写的 laser.dat 3D 图头 (就地转换, 仅头 3 行计数取整)。

    DESY Plasma_Example_2 的 laser.dat 头三行为 (n, min, spacing),
    但 n 写成浮点形式 (8.1e+01)。实测 (macOS Apple Silicon ASTRA
    构建):
      * 浮点计数 8.1e+01            -> "Error while reading file"
      * 逐值网格行 (n 后跟 n 个值)   -> 读取 3D 图时 SIGSEGV
      * 整数计数 + (n,min,spacing)  -> 正常 (本函数的输出)

    只把计数转成整数, 网格与数据逐字节不动 (物理等价)。注意: 本函数
    就地覆盖原文件 (2026-08 审计 P3), 调用前请自行备份。
    """
    path = Path(path)
    warnings.warn(
        "fix_laser_map_header 就地改写文件 %s (仅头 3 行计数取整, "
        "数据不动); 原始 DESY 头形式将被覆盖, 请确认已有备份" % path,
        UserWarning, stacklevel=2)
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

    固有歧义: 三个维度 n 均为 2 时, 逐值头每行恰 3 token, 与紧凑头
    (n, min, spacing) 完全同形, 且两种解析的数据长度都匹配, 无法从
    token 流区分; 此时固定按手册逐值头解析 (网格 = [min, spacing],
    而非紧凑展开的 [min, min+spacing]) 并发出 UserWarning。

    Returns:
        (x, y, z, F) — F 为 (nx, ny, nz) 数组, 按 F[ix, iy, iz] 索引,
        SI 单位按文件名约定 (ex/ey/ez: V/m; bx/by/bz: T; 数值原样返回)。
    """
    path = Path(path)
    # 批 5: 单次读入, 全局 token 流 (此前紧凑分支读一遍、逐值分支
    # 又 read_text() 读第二遍并逐 token 转 float, 65MB laser.dat 翻倍)
    vals = [float(v) for v in path.read_text(encoding="utf-8").split()]
    if len(vals) < 7:
        raise ValueError("3D map too short: " + str(path))

    def _per_value():
        """逐值头 (手册格式): 每维 n 后跟 n 个网格值。"""
        idx = 0
        grids = []
        for _ in range(3):
            n = int(vals[idx])
            if n < 1:
                raise ValueError("3D map grid size invalid: " + str(path))
            idx += 1
            grids.append(np.array(vals[idx:idx + n], dtype=float))
            idx += n
        return grids, np.array(vals[idx:], dtype=float)

    def _compact():
        """紧凑头 (n, min, spacing): 前 9 token 三组, 就地展开。"""
        head9 = vals[:9]
        if len(head9) < 9:
            raise ValueError("3D map compact header too short: " + str(path))
        ns = [int(head9[0]), int(head9[3]), int(head9[6])]
        if any(n < 1 for n in ns):
            raise ValueError("3D map grid size invalid: " + str(path))
        grids = [np.array([head9[1] + j * head9[2] for j in range(ns[0])],
                          dtype=float),
                 np.array([head9[4] + j * head9[5] for j in range(ns[1])],
                          dtype=float),
                 np.array([head9[7] + j * head9[8] for j in range(ns[2])],
                          dtype=float)]
        return grids, np.array(vals[9:], dtype=float)

    # 优先逐值头 (手册格式); 仅当逐值头解析失败或网格积与剩余
    # 数据长度不匹配时回退紧凑头 (MATLAB/DESY 激光图变体)。逐值头
    # 各维 n=2 时每行恰 3 token, 与紧凑头同形, 此前的"紧凑优先"
    # 会把 (2,2,2) 逐值头误判成紧凑头并静默错轴; 反向的 (2,2,2)
    # 紧凑头同样无法区分 (固有歧义), 固定按逐值头解析并显式告警。
    grids = data = None
    try:
        g, d = _per_value()
        if len(d) == len(g[0]) * len(g[1]) * len(g[2]):
            grids, data = g, d
    except (ValueError, IndexError):
        pass
    if grids is None:
        grids, data = _compact()
    elif all(len(g) == 2 for g in grids):
        warnings.warn(
            "3D 场图 %s: 各维 n=2 时逐值头与紧凑头 (n,min,spacing) "
            "token 完全同形且数据长度均匹配 (固有歧义), 已按手册逐值头"
            "解析: 网格 = [min, spacing] (紧凑展开则为 [min, min+spacing])"
            % path, UserWarning, stacklevel=2)

    x, y, z = grids
    nx, ny, nz = len(x), len(y), len(z)
    if len(data) < nx * ny * nz:
        raise ValueError("3D map data truncated: " + str(path))
    f = data[:nx * ny * nz].reshape(nx, ny, nz, order="F")
    return x, y, z, f


@dataclass
class FieldMap3D:
    """三维场图 (三分量 + 单位元数据, 矢量/等值绘图入口).

    网格与分量均为 SI: 位置 [m], 场数值按文件名约定 (bx/by/bz -> T,
    ex/ey/ez -> V/m, 数值原样返回)。缺失的分量按全零处理。
    """

    x: np.ndarray
    y: np.ndarray
    z: np.ndarray
    fx: np.ndarray
    fy: np.ndarray
    fz: np.ndarray
    unit: str = ""
    quantity: str = "F"   # 'B' / 'E' / 'F', 用于标签如 |B| [T]
    source: str = ""      # 来源路径 (告警/调试用)

    @property
    def magnitude(self) -> np.ndarray:
        """矢量模长 |F| = sqrt(fx^2 + fy^2 + fz^2)。"""
        return np.sqrt(self.fx ** 2 + self.fy ** 2 + self.fz ** 2)

    def component(self, name: str) -> np.ndarray:
        """按 'x'/'y'/'z' 取分量数组。"""
        return {"x": self.fx, "y": self.fy, "z": self.fz}[name]


# 分量后缀 -> (FieldMap3D 字段, 单位, 物理量字母)
# 电分量在前: 主名不带后缀且两族并存时 (如 Cavity_Example 的 3D_test)
# 优先按电场 (V/m) 解释; 只有磁分量时回退磁场 (T)。
_FIELD_COMPONENT_SUFFIXES = {
    ".ex": ("fx", "V/m", "E"), ".ey": ("fy", "V/m", "E"),
    ".ez": ("fz", "V/m", "E"),
    ".bx": ("fx", "T", "B"), ".by": ("fy", "T", "B"), ".bz": ("fz", "T", "B"),
}


def read_3d_field_map_components(base):
    """读取 3D 场图全部分量并组装 FieldMap3D (矢量绘图的数据入口).

    Args:
        base: 主名 (如 3D_Dipole) 或任一分量文件 (如 3D_Dipole.by);
              按文件名约定推导单位 (bx/by/bz -> T, ex/ey/ez -> V/m)。

    Returns:
        FieldMap3D — 三分量数组 (nx, ny, nz), 缺失的分量文件按全零。

    Raises:
        ValueError: 找不到任何分量文件, 或分量网格不一致。

    注意: 电/磁分量并存时只取首个出现的族 (suffix 顺序 ex..bz, 故
    电场优先), 避免磁场文件覆盖电场分量 (如 Cavity_Example/3D_test
    同时含 .ex/.ey/.ez 与 .bx/.by/.bz 六文件)。
    """
    base = Path(base)
    unit, quantity, explicit, family = "", "F", False, None
    for suffix, (_, u, q) in _FIELD_COMPONENT_SUFFIXES.items():
        if base.name.lower().endswith(suffix):
            base = base.with_name(base.name[:-len(suffix)])
            unit, quantity = u, q
            explicit = True
            family = q
            break
    grids = None
    comps = {}
    other_family_seen = False
    for suffix, (key, u, q) in _FIELD_COMPONENT_SUFFIXES.items():
        p = base.with_name(base.name + suffix)
        if not p.exists():
            continue
        if family is None:
            # 主名无后缀: 以实际找到的分量文件后缀推断单位/族
            family = q
            if not explicit:
                unit, quantity = u, q
                explicit = True
        if q != family:
            other_family_seen = True   # 电/磁并存: 只取首族 (电场优先)
            continue
        gx, gy, gz, f = read_3d_field_map(p)
        if grids is None:
            grids = (gx, gy, gz)
        elif not (len(grids[0]) == len(gx) and len(grids[1]) == len(gy)
                  and len(grids[2]) == len(gz)
                  and np.allclose(grids[0], gx) and np.allclose(grids[1], gy)
                  and np.allclose(grids[2], gz)):
            # 手册允许各分量网格不同 (仅首末 z 平面一致); 插值到
            # 首分量网格 (2026-08 审计 P3: 旧实现直接拒绝)
            warnings.warn(
                "3D 场图 %s: 分量 %s 网格与其余分量不一致 (手册允许), "
                "已线性插值到公共网格" % (base, suffix), UserWarning,
                stacklevel=2)
            from scipy.interpolate import RegularGridInterpolator
            interp = RegularGridInterpolator(
                (gx, gy, gz), f, bounds_error=False, fill_value=0.0)
            Xg, Yg, Zg = np.meshgrid(grids[0], grids[1], grids[2],
                                     indexing="ij")
            f = interp(np.column_stack([Xg.ravel(), Yg.ravel(), Zg.ravel()])
                       ).reshape(Xg.shape)
        comps[key] = f
    if grids is None:
        raise ValueError("3D 场图分量文件缺失: " + str(base))
    if other_family_seen:
        warnings.warn(
            "3D 场图 %s: 同时存在电场族与磁场族分量文件, 仅显示 %s 族 "
            "(%s 族未显示)" % (base, quantity, "B" if quantity == "E" else "E"),
            UserWarning, stacklevel=2)
    x, y, z = grids
    shape = (len(x), len(y), len(z))
    zeros = np.zeros(shape)
    return FieldMap3D(
        x=x, y=y, z=z,
        fx=comps.get("fx", zeros),
        fy=comps.get("fy", zeros),
        fz=comps.get("fz", zeros),
        unit=unit, quantity=quantity, source=str(base))

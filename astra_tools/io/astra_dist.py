"""ASTRA phase-space (particle distribution) reader.

Formats per ASTRA Manual V3.2, Table 1:

Binary:  TWO layouts are supported.
         (a) Real ASTRA output (V3.2/V4.0, verified against the local
         Apple-Silicon build): Fortran sequential unformatted records,
         one per particle: [i32 record_len=72][x, y, z, px, py, pz,
         clock, charge as 8 float64][species, status as 2 int32]
         [i32 record_len=72]. The FIRST record is the reference
         particle in ABSOLUTE coordinates; the remaining records hold
         z/pz/clock RELATIVE to it (manual Table 1, same semantics as
         the ASCII format). 500-particle file = 40,000 bytes.
         (b) Legacy stream layout written by write_distribution /
         constructed test files: 5x float64 header (ref time [ns],
         ref momentum [eV/c], total charge [nC], x_ref [m], y_ref [m])
         followed by N particles of 10 float64 (manual Table 1 is
         always 10 columns):
             x, y, z, px, py, pz, clock, charge, species, status
         Units: x/y/z [m], px/py/pz [eV/c], clock [ns], charge [nC].
         species (col 9) = particle species index (1=e-, 2=e+, 3=p,
         4=H+, 5..14=user defined; manual Table 1) - NOT a running
         particle number. Legacy 9-column files (x..status, no species
         column) are still readable with a warning.

ASCII:   no header line; the first particle row is the reference
         particle in ABSOLUTE coordinates, the remaining rows have
         z/pz/clock RELATIVE to it (Fortran format 1P,8E12.4,2I4;
         High_res writes 1P,8E20.12,2I4). Legacy project files with a
         5-value header line are also accepted.

The 9/10-column ambiguity (N divisible by both 9 and 10) is resolved
by column semantics, not by any running-index heuristic: a 10-column
reading is accepted when the species column holds integers in [1..14]
and the status column holds small integers.

Units on output are canonical SI (clock converted ns -> s).
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional

import numpy as np

from ..constants import NS_TO_S
from ..distribution import Distribution

logger = logging.getLogger(__name__)

# File suffixes that are NOT phase-space files (ASTRA Manual Table 3/4)
_NON_PHASE_SUFFIXES = {
    ".xemit", ".yemit", ".zemit", ".cemit", ".c99emit", ".tremit",
    ".xemit2", ".yemit2", ".cr_emit", ".sub_emit", ".sigma", ".ref",
    ".log", ".track", ".cathode", ".larmor", ".density", ".landf",
    ".pscan", ".scan", ".error", ".tstep", ".fields", ".tcheck",
}

_PHASE_SUFFIXES = {".ini", ".ast", ".inp", ".zpos"}


class AstraDistributionReader:
    """Reader for ASTRA binary / ASCII particle distribution files."""

    format_name = "astra"

    def probe(self, path: Path) -> bool:
        """Heuristic format detection.

        Rejects emit/sigma/ref/log files by suffix, accepts .ini/.ast/.inp
        and numeric suffixes (.001), and for other files (e.g. .dat, which
        may be a field map) requires the content to actually look like a
        particle file.
        """
        path = Path(path)
        suffixes_lower = {s.lower() for s in path.suffixes}
        if suffixes_lower & _NON_PHASE_SUFFIXES:
            return False
        if path.suffix.lower() in _PHASE_SUFFIXES:
            return True
        # Numeric suffix like .001 or .0050.001 (ASTRA tracking output)
        if path.suffix.lstrip(".").isdigit():
            return True
        # Anything else (e.g. .dat): only accept if content looks like
        # particles. Field maps (2-column tables) must be rejected.
        # 先按字节判断文本/二进制: ASCII 场图的前几个字节读成 float64
        # 会变成 denormal 小浮点并通过二进制头的范围检查 (误判)
        with open(path, "rb") as f:
            head_bytes = f.read(4096)
        if b"\x00" not in head_bytes:
            # 文本: 只走 ASCII 探测
            ncols = self._probe_ascii_ncols(path)
            if ncols is not None:
                if ncols >= 9:
                    return True
                if ncols == 5 and self._probe_ascii_second_line(path) in (9, 10):
                    return True
            # 内容不像 ASCII 数值表: 小二进制文件可能全无 0x00 字节
            # (2026-08 审计 P3), 尝试二进制探测再放弃
            try:
                self._probe_binary(path)
                return True
            except Exception:
                return False
        try:
            self._probe_binary(path)
            return True
        except Exception:
            pass
        try:
            ncols = self._probe_ascii_ncols(path)
            if ncols is None:
                return False
            if ncols >= 9:
                return True
            # 5 值头行 + 9/10 列粒子行的 ASCII 文件 (第二行才是粒子行)
            return ncols == 5 and self._probe_ascii_second_line(path) in (9, 10)
        except Exception:
            return False

    @staticmethod
    def _probe_ascii_second_line(path: Path):
        """ASCII 文件第二行 (跳过注释/空行) 的列数; 无则 None。"""
        with open(path, "rb") as f:
            head = f.read(8192)
        text = head.decode("ascii", errors="ignore")
        rows = [ln.split() for ln in text.splitlines()
                if ln.strip() and not ln.strip().startswith(("#", "!"))]
        if len(rows) < 2:
            return None
        return len(rows[1])

    # -- low-level probes ---------------------------------------------

    @staticmethod
    def _probe_binary(path: Path) -> None:
        data = np.fromfile(path, dtype=np.float64)
        if len(data) < 5:
            raise ValueError("file too small for ASTRA header")
        if not AstraDistributionReader._binary_plausible(data):
            swapped = np.ascontiguousarray(data.byteswap())
            if not AstraDistributionReader._binary_plausible(swapped):
                # 真实 ASTRA 二进制 = Fortran unformatted 记录流
                # (非 5 值头流式布局), 记录探测通过才算二进制。
                if not (AstraDistributionReader._records_plausible(path)
                        or AstraDistributionReader._records_plausible(
                            path, swap=True)):
                    raise ValueError(
                        "binary header/body out of physical range")

    @staticmethod
    def _binary_plausible(data: np.ndarray) -> bool:
        """二进制解释合理性: 头 5 值在物理范围内, 且 body 大部分是
        正常浮点 (字节序错误时几乎全是 denormal/inf/nan)."""
        if len(data) < 5:
            return False
        hdr = data[:5]
        if not (np.all(np.isfinite(hdr)) and abs(hdr[0]) <= 1e6
                and -1.0 <= hdr[1] <= 1e13 and abs(hdr[2]) <= 1e6):
            return False
        body = data[5:]
        if len(body) == 0:
            return True
        with np.errstate(invalid="ignore", over="ignore"):
            normal = np.isfinite(body) & (np.abs(body) >= 1e-300)
        return float(np.mean(normal)) >= 0.5

    @staticmethod
    def _probe_ascii_ncols(path: Path):
        with open(path, "rb") as f:
            head = f.read(4096)
        if b"\x00" in head:
            return None  # binary file
        text = head.decode("ascii", errors="ignore")
        first = text.splitlines()[0] if text.splitlines() else ""
        parts = first.split()
        return len(parts) if parts else None

    # -- public API ---------------------------------------------------

    def read(self, path: Path) -> Distribution:
        """Read an ASTRA distribution file (binary or ASCII, auto-detected)."""
        path = Path(path)
        with open(path, "rb") as f:
            head = f.read(4096)
        if b"\x00" in head:
            return self._read_binary(path)
        # 无 0x00 字节: 内容必须真的像 ASCII 数值表才走 ASCII 路径;
        # 小二进制文件可能全无 0x00 字节 (2026-08 审计 P3)。
        try:
            ncols = self._probe_ascii_ncols(path)
            if ncols in (5, 9, 10):
                return self._read_ascii(path)
        except Exception:
            pass
        return self._read_binary(path)

    # -- binary -------------------------------------------------------

    # 真实 ASTRA 二进制相位空间文件 (V3.2/V4.0, macOS gfortran 实测):
    # Fortran sequential unformatted 记录流, 每条记录
    #   [i32 len=72][8×f64: x y z px py pz clock q][i32 species][i32 status][i32 len=72]
    # 共 80 字节; 首条记录 = 参考粒子绝对坐标, 其余记录 z/pz/clock 相对。
    # (R1 golden: examples/Manual_Example/golden/Example_binary.001,
    #  500 粒子 = 40000 字节; 与同 run ASCII 相位 dump 逐粒子一致。)
    _RECORD_SIZE = 80
    _RECORD_LEN = 72

    @staticmethod
    def _records_plausible(path: Path, swap: bool = False) -> bool:
        """记录流布局探测: 80 字节记录 + 长度标记 + 载荷物理范围."""
        try:
            raw = Path(path).read_bytes()
        except OSError:
            return False
        n = len(raw)
        if n == 0 or n % AstraDistributionReader._RECORD_SIZE != 0:
            return False
        dt = np.dtype([("mlen", "<i4"), ("vals", "<8f8"), ("idx", "<i4"),
                       ("stat", "<i4"), ("tlen", "<i4")])
        if swap:
            dt = dt.newbyteorder("S")
        rec = np.frombuffer(raw, dtype=dt)
        mlen = np.asarray(rec["mlen"])
        tlen = np.asarray(rec["tlen"])
        if not (bool(np.all(mlen == AstraDistributionReader._RECORD_LEN))
                and bool(np.all(tlen == AstraDistributionReader._RECORD_LEN))):
            return False
        p = np.asarray(rec["vals"])
        idx = np.asarray(rec["idx"])
        stat = np.asarray(rec["stat"])
        with np.errstate(invalid="ignore", over="ignore"):
            ok_finite = bool(np.all(np.isfinite(p)))
            ok_xy = bool(np.all(np.abs(p[:, :2]) <= 1e3))
            ok_z = bool(np.all(np.abs(p[:, 2]) <= 1e4))
            ok_p = bool(np.all(np.abs(p[:, 3:6]) <= 1e15))
            ok_cq = bool(np.all(np.abs(p[:, 6:8]) <= 1e6))
        return (ok_finite and ok_xy and ok_z and ok_p and ok_cq
                and bool(np.all((idx >= 1) & (idx <= 14)))
                and bool(np.all(np.abs(stat) <= 1e6)))

    def _read_records(self, path: Path, swap: bool = False):
        """解析真实 ASTRA unformatted 记录流; 非记录流返回 None."""
        if not self._records_plausible(path, swap=swap):
            return None
        raw = Path(path).read_bytes()
        dt = np.dtype([("mlen", "<i4"), ("vals", "<8f8"), ("idx", "<i4"),
                       ("stat", "<i4"), ("tlen", "<i4")])
        if swap:
            dt = dt.newbyteorder("S")
        rec = np.frombuffer(raw, dtype=dt)
        p = np.asarray(rec["vals"])        # (N, 8): x y z px py pz clock q
        idx = np.asarray(rec["idx"]).astype(np.int32)
        stat = np.asarray(rec["stat"]).astype(np.int32)

        # 首条记录 = 参考粒子绝对坐标; 其余记录 z/pz/clock 相对 (Table 1)
        ref = p[0]
        z_abs = p[:, 2].copy()
        z_abs[1:] += ref[2]
        pz_abs = p[:, 5].copy()
        pz_abs[1:] += ref[5]
        clock_abs = p[:, 6].copy()
        clock_abs[1:] += ref[6]

        dist = Distribution(
            x=p[:, 0], y=p[:, 1], z=z_abs,
            px=p[:, 3], py=p[:, 4], pz=pz_abs,
            clock=clock_abs * NS_TO_S,      # [s]
            charge=p[:, 7],                  # nC
            status=stat,
            index=idx,
            ref_time_ns=float(ref[6]),
            ref_momentum_eVc=float(ref[5]),
            total_charge_nC=float(np.abs(np.sum(p[:, 7]))),   # |Q| 约定
            ref_x_m=float(ref[0]),
            ref_y_m=float(ref[1]),
            ref_z_m=float(ref[2]),
            source=str(path),
            format="astra_binary_records",
        )
        logger.info(
            "Read ASTRA binary (record stream) %s: p_ref=%.3f MeV/c, "
            "Q=%.4f nC, N=%d (%d active)",
            path.name, ref[5] * 1e-6, dist.active_charge_nC,
            dist.n_particle, dist.n_active,
        )
        return dist

    def _read_binary(self, path: Path) -> Distribution:
        # 真实 ASTRA unformatted 记录流优先 (2026-08 R1 真跑发现:
        # 官方二进制输出并非 5 值头流式布局)。
        rec = self._read_records(path)
        if rec is None and AstraDistributionReader._records_plausible(
                path, swap=True):
            rec = self._read_records(path, swap=True)
        if rec is not None:
            return rec
        data = np.fromfile(path, dtype=np.float64)
        if len(data) < 5:
            raise ValueError(
                "file " + str(path) + " too small for ASTRA header: "
                + str(len(data)) + " values"
            )

        # 字节序探测 (2026-08 审计 P3): 原生序不合理时试大端。
        if not self._binary_plausible(data):
            swapped = np.ascontiguousarray(data.byteswap())
            if self._binary_plausible(swapped):
                data = swapped
                logger.warning("binary file %s read with swapped byte order", path.name)
            else:
                raise ValueError(
                    "binary file " + str(path.name) + " has implausible header/body "
                    "(not an ASTRA distribution, or corrupted)"
                )

        header = data[:5].copy()
        body = data[5:]

        n9, r9 = divmod(len(body), 9)
        n10, r10 = divmod(len(body), 10)

        # 手册 Table 1: 二进制恒为 10 列 (含 species)。9 列是遗留非标准
        # 格式; 两者长度均可整除时按列语义判定 (species/status), 不用
        # 任何"递增编号"启发式 (2026-08 审计 P1-1: 第 9 列是粒子种类,
        # 电子束恒为 1)。
        if r10 == 0:
            test = body[: n10 * 10].reshape(n10, 10)
            if self._species_status_plausible(test):
                n_particles, n_cols = n10, 10
            elif r9 == 0:
                n_particles, n_cols = n9, 9
                logger.warning(
                    "binary file %s read as legacy 9-column format "
                    "(no species column; ASTRA writes 10 columns)", path.name)
            else:
                raise ValueError(
                    "file " + str(path.name) + " has " + str(len(body))
                    + " body values divisible by 10 but the species/status "
                    "columns are implausible"
                )
        elif r9 == 0:
            n_particles, n_cols = n9, 9
            logger.warning(
                "binary file %s read as legacy 9-column format "
                "(no species column; ASTRA writes 10 columns)", path.name)
        else:
            raise ValueError(
                "file " + str(path.name) + " has ambiguous size: "
                + str(len(body)) + " values is neither 9 nor 10 columns "
                "(remainders " + str(r9) + "/" + str(r10) + ")"
            )

        p = body[: n_particles * n_cols].reshape(n_particles, n_cols)

        # Column layout: x y z px py pz clock charge [species] status
        # Manual Table 1: z, pz and clock are RELATIVE to the reference
        # particle -> convert to absolute using the header values.
        status_col = 9 if n_cols == 10 else 8
        index = p[:, 8].astype(np.int32) if n_cols == 10 else None
        status_f = p[:, status_col]
        if np.any(np.abs(status_f) > 2.147e9):
            raise ValueError(
                "status values out of int32 range in " + str(path.name)
                + " (corrupt file?)"
            )
        pz_abs = p[:, 5] + float(header[1])
        clock_abs = (p[:, 6] + float(header[0])) * NS_TO_S
        # z offset of the reference particle is not stored in the binary
        # header: keep z relative (ref_z_m stays 0).
        dist = Distribution(
            x=p[:, 0], y=p[:, 1], z=p[:, 2],
            px=p[:, 3], py=p[:, 4], pz=pz_abs,
            clock=clock_abs,          # [s]
            charge=p[:, 7],            # nC
            status=status_f.astype(np.int32),
            index=index,
            ref_time_ns=float(header[0]),
            ref_momentum_eVc=float(header[1]),
            total_charge_nC=float(header[2]),
            ref_x_m=float(header[3]) if len(header) > 3 else 0.0,
            ref_y_m=float(header[4]) if len(header) > 4 else 0.0,
            source=str(path),
            format="astra_binary",
        )

        logger.info(
            "Read ASTRA binary %s: p_ref=%.3f MeV/c, Q=%.4f nC, "
            "N=%d (%d active)",
            path.name, header[1] * 1e-6, dist.active_charge_nC,
            dist.n_particle, dist.n_active,
        )
        return dist

    # -- ASCII --------------------------------------------------------

    def _read_ascii(self, path: Path) -> Distribution:
        # 批 5: 向量化 (np.loadtxt), 不再逐行 split+float 循环;
        # 5 值头行与 9/10 列粒子行不齐, 先探测首内容行决定 skiprows。
        # 复核修正: 跳过前导空行与 #/! 注释行 (原实现先过滤, 不可回归)。
        raw_skip = 0
        first_parts = None
        with open(path, encoding="utf-8") as f:
            for ln in f:
                s = ln.strip()
                if not s or s.startswith(("#", "!")):
                    raw_skip += 1
                    continue
                first_parts = s.split()
                break
        if first_parts is None:
            raise ValueError("empty ASTRA ASCII file: " + str(path))
        # A header line has exactly 5 values; a particle row has 9 or 10.
        has_header = len(first_parts) == 5
        data_start = 1 if has_header else 0

        header = np.zeros(5)
        if has_header:
            header[: min(5, len(first_parts))] = [float(v) for v in first_parts[:5]]

        try:
            p = np.loadtxt(path, ndmin=2, skiprows=raw_skip + data_start,
                           comments=("#", "!"), encoding="utf-8")
        except ValueError as e:
            raise ValueError(
                "malformed particle row in " + str(path.name)
                + ": expected 9-10 columns (" + str(e) + ")")
        n = p.shape[0]
        ncols = p.shape[1]
        if ncols not in (9, 10):
            raise ValueError(
                "malformed particle row in " + str(path.name)
                + ": expected 9-10 columns, got " + str(ncols))
        status_col = 9 if ncols == 10 else 8
        index = p[:, 8].astype(np.int32) if ncols == 10 else None

        # No header line: the first particle row IS the reference
        # particle with ABSOLUTE coordinates (ASTRA Manual V3.2:
        # "The first line of the file defines the coordinates of the
        # reference particle in absolute coordinates").
        ref_z_m = 0.0
        if not has_header and n > 0:
            header[0] = p[0, 6]      # clock [ns] -> ref time
            header[1] = p[0, 5]      # pz [eV/c] -> ref momentum
            header[2] = float(np.sum(p[:, 7]))  # total charge [nC]
            header[3] = p[0, 0]
            header[4] = p[0, 1]
            ref_z_m = float(p[0, 2])
            # Remaining rows: z, pz, clock are relative to the reference
            # particle (Table 1) -> shift to absolute.
            p[1:, 2] += ref_z_m
            p[1:, 5] += float(header[1])
            p[1:, 6] += float(header[0])

        # With an explicit 5-value header line the particle rows are the
        # bunch; pz/clock are still relative to the reference values
        # (z offset unknown, stays relative).
        if has_header:
            p[:, 5] += float(header[1])
            p[:, 6] += float(header[0])

        dist = Distribution(
            x=p[:, 0], y=p[:, 1], z=p[:, 2],
            px=p[:, 3], py=p[:, 4], pz=p[:, 5],
            clock=p[:, 6] * NS_TO_S,   # ns -> s
            charge=p[:, 7],            # nC
            status=p[:, status_col].astype(np.int32),
            index=index,
            ref_time_ns=float(header[0]),
            ref_momentum_eVc=float(header[1]),
            total_charge_nC=float(header[2]),
            ref_x_m=float(header[3]),
            ref_y_m=float(header[4]),
            ref_z_m=ref_z_m,
            source=str(path),
            format="astra_ascii",
        )

        logger.info(
            "Read ASTRA ASCII %s: N=%d (%d active)",
            path.name, dist.n_particle, dist.n_active,
        )
        return dist

    @staticmethod
    def _species_status_plausible(p: np.ndarray) -> bool:
        """10 列解释的列语义判定 (手册 Table 1).

        第 9 列 (species) 必须是 [1..14] 内的整数 (1=电子, 2=正电子,
        3=质子, 4=氢离子, 5-14=用户定义); 第 10 列 (status) 必须是
        |v| <= 1e6 的整数。9 列遗留文件被误排成 10 列时, species 位置
        是下一行粒子的 x (非整数), status 位置是 y (非整数) — 判定失败。
        """
        if p.shape[0] == 0:
            return True
        species = p[:, 8]
        status = p[:, 9]
        with np.errstate(invalid="ignore"):
            spec_ok = (
                bool(np.all(np.isfinite(species)))
                and bool(np.all(species == np.round(species)))
                and bool(np.all((species >= 1) & (species <= 14)))
            )
            stat_ok = (
                bool(np.all(np.isfinite(status)))
                and bool(np.all(status == np.round(status)))
                and bool(np.all(np.abs(status) <= 1e6))
            )
        return spec_ok and stat_ok


def write_distribution(
    dist: Distribution,
    path,
    format: str = "binary",
    include_index: bool = True,
    ref_z_m: Optional[float] = None,
) -> str:
    """写 ASTRA 分布文件 (postpro 5.6.4 保存新分布供继续追踪).

    与 AstraDistributionReader 的读取约定互逆 (手册 Table 1):
      * binary: 5-float64 头 (ref_time_ns, ref_momentum_eVc, |Q|,
        ref_x_m, ref_y_m) + 每粒子 10 列
        x, y, z_rel, px, py, pz_rel, clock_rel_ns, charge_nC,
        species, status。头 Q 取 |Q| (手册: 电荷符号不相关)。
      * ascii: 无头行, 第一行是参考粒子 (粒子 0) 绝对坐标 (10 列),
        其余行 z/pz/clock 相对参考粒子 — 真实 ASTRA 格式。

    Args:
        dist: Distribution。
        path: 输出路径。
        format: 'binary' (默认, ASTRA 标准输入) 或 'ascii'。
        include_index: 是否写种类列 (手册 Table 1 第 9 列; 恒写 10 列
            为 ASTRA 标准; 设 False 得到 9 列遗留格式, 仅供兼容)。
            index 为 None 时种类默认写 1 (电子), 不伪造 1..N 编号。
        ref_z_m: 二进制路径参考粒子绝对 z [m]; 默认 dist.ref_z_m
            (二进制读取时其为 0, z 即相对值, 原样写回)。ASCII 路径
            以粒子 0 为参考粒子 (真实 ASTRA 语义), 不使用此参数。

    Returns:
        写入的文件路径字符串。
    """
    path = Path(path)
    if ref_z_m is None:
        ref_z_m = dist.ref_z_m
    n = dist.n_particle
    if n == 0:
        raise ValueError("cannot write an empty distribution")
    p_ref = dist.ref_momentum_eVc or dist.mean_pz_eVc
    ref_time = dist.ref_time_ns
    z_rel = dist.z - ref_z_m
    pz_rel = dist.pz - p_ref
    clock_rel_ns = dist.clock * (1.0 / NS_TO_S) - ref_time

    if format == "ascii":
        # 真实 ASTRA ASCII 格式: 无头行; 第一行 = 参考粒子 (粒子 0)
        # 绝对坐标, 其余行 z/pz/clock 相对参考粒子 (手册 Table 1)。
        # (2026-08 审计 P2-1: 旧的 5 值头形式 ASTRA 无法解析。)
        idx0 = dist.index[0] if dist.index is not None else 1
        ref_row = [dist.x[0], dist.y[0], dist.z[0],
                   dist.px[0], dist.py[0], dist.pz[0],
                   dist.clock[0] / NS_TO_S, dist.charge[0],
                   idx0, dist.status[0]]
        with open(path, "w") as fh:
            fh.write(" ".join("%.12g" % v for v in ref_row) + "\n")
            for i in range(1, n):
                idx = dist.index[i] if dist.index is not None else 1
                row = [dist.x[i], dist.y[i], dist.z[i] - dist.z[0],
                       dist.px[i], dist.py[i], dist.pz[i] - dist.pz[0],
                       dist.clock[i] / NS_TO_S - dist.clock[0] / NS_TO_S,
                       dist.charge[i], idx, dist.status[i]]
                fh.write(" ".join("%.12g" % v for v in row) + "\n")
    else:
        # 头 Q 取 |Q| (2026-08 审计 P1-2: ASTRA 约定头电荷为正)。
        q_abs = float(np.sum(np.abs(dist.charge)))
        header = np.array([ref_time, p_ref, q_abs,
                           dist.ref_x_m, dist.ref_y_m], dtype=np.float64)
        cols = [dist.x, dist.y, z_rel, dist.px, dist.py, pz_rel,
                clock_rel_ns, dist.charge]
        if include_index:
            idx = dist.index if dist.index is not None else np.ones(n, dtype=np.int32)
            cols.append(idx.astype(np.float64))
        cols.append(dist.status.astype(np.float64))
        body = np.column_stack(cols).astype(np.float64).ravel()
        with open(path, "wb") as fh:
            header.tofile(fh)
            body.tofile(fh)
    logger.info("Wrote distribution %s (%s, N=%d)", path, format, n)
    return str(path)

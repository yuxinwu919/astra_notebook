"""ASTRA phase-space (particle distribution) reader.

Formats per ASTRA Manual V3.2, Table 1:

Binary:  5x float64 header (ref time [ns], ref momentum [eV/c], total
         charge [nC], x_ref [m], y_ref [m]) followed by N particles of
         9 or 10 float64:
             x, y, z, px, py, pz, clock, charge, [index,] status
         Units: x/y/z [m], px/py/pz [eV/c], clock [ns], charge [nC].

ASCII:   optional 5-value header line, then N rows of 10 columns
         (Fortran format 1P,8E12.4,2I4; High_res writes 1P,8E20.12,2I4).
         Generator output has NO header line; the first particle row is
         the reference particle (absolute coordinates).

The 9/10-column ambiguity in binary files is resolved by inspecting the
10th column: if it looks like sequential particle indices it is kept,
otherwise the file is read as 9 columns.

Units on output are canonical SI (clock converted ns -> s).
"""

from __future__ import annotations

import logging
from pathlib import Path

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
            if ncols is None:
                return False
            if ncols >= 9:
                return True
            return ncols == 5 and self._probe_ascii_second_line(path) in (9, 10)
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
        data = np.fromfile(path, dtype=np.float64, count=5)
        if len(data) < 5:
            raise ValueError("file too small for ASTRA header")
        if abs(data[0]) > 1e6 or data[1] < -1 or data[1] > 1e13:
            raise ValueError("header values out of physical range")

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
        is_ascii = b"\x00" not in head
        if is_ascii:
            return self._read_ascii(path)
        return self._read_binary(path)

    # -- binary -------------------------------------------------------

    def _read_binary(self, path: Path) -> Distribution:
        data = np.fromfile(path, dtype=np.float64)
        if len(data) < 5:
            raise ValueError(
                "file " + str(path) + " too small for ASTRA header: "
                + str(len(data)) + " values"
            )

        header = data[:5].copy()
        body = data[5:]

        n9, r9 = divmod(len(body), 9)
        n10, r10 = divmod(len(body), 10)

        if r9 == 0 and r10 == 0:
            # Ambiguous: inspect column 8 (particle index; column 9 is
            # status and is ~constant, so it can never look sequential)
            test = body[: n10 * 10].reshape(n10, 10)
            if self._looks_like_index(test[:, 8]):
                n_particles, n_cols = n10, 10
            else:
                n_particles, n_cols = n9, 9
        elif r9 == 0:
            n_particles, n_cols = n9, 9
        elif r10 == 0:
            n_particles, n_cols = n10, 10
            logger.warning("10-column format detected for %s (with particle index)", path.name)
        else:
            raise ValueError(
                "file " + str(path.name) + " has ambiguous size: "
                + str(len(body)) + " values is neither 9 nor 10 columns "
                "(remainders " + str(r9) + "/" + str(r10) + ")"
            )

        p = body[: n_particles * n_cols].reshape(n_particles, n_cols)

        # Column layout: x y z px py pz clock charge [index] status
        # Manual Table 1: z, pz and clock are RELATIVE to the reference
        # particle -> convert to absolute using the header values.
        status_col = 9 if n_cols == 10 else 8
        index = p[:, 8].astype(np.int32) if n_cols == 10 else None
        pz_abs = p[:, 5] + float(header[1])
        clock_abs = (p[:, 6] + float(header[0])) * NS_TO_S
        # z offset of the reference particle is not stored in the binary
        # header: keep z relative (ref_z_m stays 0).
        dist = Distribution(
            x=p[:, 0], y=p[:, 1], z=p[:, 2],
            px=p[:, 3], py=p[:, 4], pz=pz_abs,
            clock=clock_abs,          # [s]
            charge=p[:, 7],            # nC
            status=p[:, status_col].astype(np.int32),
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
        with open(path, encoding="utf-8") as f:
            lines = [ln.strip() for ln in f if ln.strip() and not ln.startswith(("#", "!"))]
        if not lines:
            raise ValueError("empty ASTRA ASCII file: " + str(path))

        first_parts = lines[0].split()
        # A header line has exactly 5 values; a particle row has 9 or 10.
        has_header = len(first_parts) == 5
        data_start = 1 if has_header else 0

        header = np.zeros(5)
        if has_header:
            header[: min(5, len(first_parts))] = [float(v) for v in first_parts[:5]]

        rows = []
        for ln in lines[data_start:]:
            vals = ln.split()
            if len(vals) not in (9, 10):
                raise ValueError(
                    "malformed particle row in " + str(path.name)
                    + ": expected 9-10 columns, got " + str(len(vals))
                )
            rows.append([float(v) for v in vals])

        p = np.asarray(rows, dtype=float)
        n = p.shape[0]
        ncols = p.shape[1]
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
    def _looks_like_index(col: np.ndarray) -> bool:
        """True if the column looks like sequential particle indices."""
        if len(col) < 2:
            return False
        diffs = np.diff(col)
        return bool(np.all(np.abs(diffs - 1.0) < 0.1))

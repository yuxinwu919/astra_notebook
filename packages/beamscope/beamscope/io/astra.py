"""ASTRA format reader.

Handles ASTRA binary and ASCII particle distribution files.
Reference: ASTRA Manual V3.2, Table 1.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

from ..constants import NS_TO_S
from ..distribution import Distribution

logger = logging.getLogger(__name__)


class AstraReader:
    """Reader for ASTRA particle distribution files.

    Supports both binary (standard) and ASCII formats.
    Binary format: 5×float64 header + N×9×float64 particles.
    ASCII format: space-separated text, optional particle_index column.
    """

    format_name = "astra"

    # ASTRA binary dtype: 9 columns (standard)
    ASTRA_DTYPE_9COL = np.dtype(
        [
            ("x", "f8"),
            ("y", "f8"),
            ("z", "f8"),
            ("px", "f8"),
            ("py", "f8"),
            ("pz", "f8"),
            ("clock", "f8"),
            ("macro_charge", "f8"),
            ("status_flag", "f8"),
        ]
    )

    # ASTRA binary dtype: 10 columns (with particle_index)
    ASTRA_DTYPE_10COL = np.dtype(
        [
            ("x", "f8"),
            ("y", "f8"),
            ("z", "f8"),
            ("px", "f8"),
            ("py", "f8"),
            ("pz", "f8"),
            ("clock", "f8"),
            ("macro_charge", "f8"),
            ("status_flag", "f8"),
            ("particle_index", "f8"),
        ]
    )

    # File name patterns that indicate emit/sigma/ref/log files (NOT phase space)
    _EMIT_SIGMA_EXTENSIONS = {'.xemit', '.yemit', '.zemit', '.cemit',
                               '.sigma', '.ref', '.log'}

    def probe(self, path: Path) -> bool:
        """Check if a file appears to be ASTRA phase space format.

        Heuristic: check file extension and attempt to read header.
        Recognizes .ini, .dat, .ast, and numeric suffixes like .001
        (ASTRA phase space output).  Rejects emit/sigma/ref/log files.
        """
        # Check intermediate extensions to reject emit/sigma/ref/log files
        suffixes_lower = {s.lower() for s in path.suffixes}
        if suffixes_lower & self._EMIT_SIGMA_EXTENSIONS:
            return False
        # Known ASTRA phase space extensions
        if path.suffix.lower() in (".ini", ".dat", ".ast", ".inp"):
            return True
        # Numeric suffix like .001, .0050.001 (ASTRA tracking output)
        if path.suffix.lstrip(".").isdigit():
            return True
        # Try reading as binary — if it looks like valid ASTRA data, accept it
        try:
            self._probe_binary(path)
            return True
        except Exception:
            pass
        return False

    def _probe_binary(self, path: Path) -> None:
        """Quick probe: check if file starts with reasonable ASTRA header values."""
        data = np.fromfile(path, dtype=np.float64, count=5)
        if len(data) < 5:
            raise ValueError("File too small")
        # Header values should be within reasonable physical ranges
        # ref_time [ns]: -1000 to 1000
        # ref_energy [eV]: 0 to 1e12
        # Q_total [nC]: -1e6 to 1e6
        if abs(data[0]) > 1e6 or data[1] < -1 or data[1] > 1e13:
            raise ValueError("Header values out of physical range")

    def read(self, path: Path) -> Distribution:
        """Read an ASTRA distribution file.

        Args:
            path: Path to .ini or .dat file.

        Returns:
            Distribution instance.
        """
        # Try ocelot first (handles all ASTRA format variants correctly)
        try:
            import ocelot.adaptors.astra2ocelot as a2o
            p_array = a2o.astraBeam2particleArray(str(path))
            return self._from_ocelot(p_array, path)
        except ImportError:
            pass
        except Exception as e:
            logger.debug(f"ocelot failed for '{path.name}': {e}, trying built-in parser...")

        # Try lume-astra second
        try:
            import lume_astra
            data_dict = lume_astra.read_astra(str(path))
            return self._from_lume(data_dict, path)
        except ImportError:
            pass
        except Exception as e:
            logger.debug(f"lume-astra failed for '{path.name}': {e}, trying built-in parser...")

        # Fall back to built-in parser
        try:
            with open(path, "rb") as f:
                head = f.read(128)
            is_ascii = b"\x00" not in head and b"E" in head
        except Exception:
            is_ascii = False

        if is_ascii:
            return self._read_ascii(path)
        return self._read_binary(path)

    def _from_ocelot(self, p_array, path: Path) -> Distribution:
        """Convert an ocelot ParticleArray to a beamscope Distribution.

        ocelot coordinate conventions:
          - .x(), .y(): transverse position [m]
          - .tau(): longitudinal position z - ct [m]
          - .px(), .py(): transverse angle x'=px/p0, y'=py/p0 [rad] (dimensionless)
          - .p(): relative momentum deviation δ = (p-p0)/p0
          - .E: reference energy [GeV]
          - .q_array: macro-particle charge [C]
        """
        n = p_array.n
        p_ref = p_array.E * 1e9  # GeV → eV/c

        delta = p_array.p()
        pz = p_ref * (1.0 + delta)
        px = p_array.px() * pz
        py = p_array.py() * pz
        q_nC = p_array.q_array * 1e9  # C → nC

        dist = Distribution(
            x=p_array.x(),
            y=p_array.y(),
            z=p_array.tau(),
            px=px,
            py=py,
            pz=pz,
            clock=np.zeros(n),
            charge=q_nC,
            status=np.zeros(n, dtype=np.int32),
            ref_momentum_eVc=p_ref,
            total_charge_nC=float(np.sum(q_nC)),
            source=str(path),
            format="astra_binary",
        )

        logger.info(
            f"Read '{path.name}' via ocelot: "
            f"p_ref={p_ref * 1e-6:.2f} MeV/c, "
            f"N={n}"
        )

        return dist

    def _from_lume(self, data_dict: dict, path: Path) -> Distribution:
        """Convert lume-astra output dict to Distribution."""
        key_map = {
            "x": "x", "y": "y", "z": "z",
            "px": "px", "py": "py", "pz": "pz",
            "t": "clock", "clock": "clock",
            "macro_charge": "charge", "charge": "charge",
            "status": "status_flag", "flag": "status_flag",
        }
        mapped: dict[str, np.ndarray] = {}
        for src, dst in key_map.items():
            if src in data_dict:
                mapped[dst] = np.asarray(data_dict[src])

        n = len(mapped.get("x", []))
        if "clock" in mapped:
            mapped["clock"] = mapped["clock"] * NS_TO_S
        else:
            mapped["clock"] = np.zeros(n)
        if "status_flag" not in mapped:
            mapped["status_flag"] = np.zeros(n)
        if "charge" not in mapped:
            mapped["charge"] = np.zeros(n)

        dist = Distribution(
            x=mapped.get("x", np.zeros(n)),
            y=mapped.get("y", np.zeros(n)),
            z=mapped.get("z", np.zeros(n)),
            px=mapped.get("px", np.zeros(n)),
            py=mapped.get("py", np.zeros(n)),
            pz=mapped.get("pz", np.zeros(n)),
            clock=mapped["clock"],
            charge=mapped["charge"],
            status=mapped["status_flag"].astype(np.int32),
            ref_momentum_eVc=float(data_dict.get("ref_energy", data_dict.get("p_ref", 0.0))),
            source=str(path),
            format="astra_binary",
        )

        logger.info(f"Read '{path.name}' via lume-astra: N={n}")
        return dist

    def _read_binary(self, path: Path) -> Distribution:
        """Read binary ASTRA format (standard)."""
        data = np.fromfile(path, dtype=np.float64)

        if len(data) < 5:
            raise ValueError(f"File '{path}' too small for ASTRA header: {len(data)} values")

        header = data[:5].copy()

        # Detect 9-col vs 10-col format
        n_particles_9 = (len(data) - 5) // 9
        n_particles_10 = (len(data) - 5) // 10
        remainder_9 = (len(data) - 5) % 9
        remainder_10 = (len(data) - 5) % 10

        if remainder_9 == 0 and remainder_10 == 0:
            # Ambiguous: prefer 9-column, inspect 10th column
            test_particles = data[5 : 5 + n_particles_10 * 10].reshape(n_particles_10, 10)
            tenth_col = test_particles[:, 9]
            # If 10th column looks like particle_index (sequential integers), it's 10-col
            if self._looks_like_index(tenth_col):
                n_particles = n_particles_10
                n_cols = 10
            else:
                n_particles = n_particles_9
                n_cols = 9
        elif remainder_9 == 0:
            n_particles = n_particles_9
            n_cols = 9
        elif remainder_10 == 0:
            n_particles = n_particles_10
            n_cols = 10
            logger.warning(f"10-column format detected for '{path.name}' (with particle_index)")
        else:
            # Best effort: try 9-column
            logger.warning(
                f"File '{path.name}' has ambiguous size: {len(data)} values "
                f"(9col remainder={remainder_9}, 10col remainder={remainder_10}). "
                f"Assuming 9-column with {n_particles_9} particles."
            )
            n_particles = n_particles_9
            n_cols = 9

        particles = data[5 : 5 + n_particles * n_cols].reshape(n_particles, n_cols)

        if n_cols == 10:
            # Strip particle_index column
            particles = particles[:, :9]

        dist = Distribution(
            x=particles[:, 0],
            y=particles[:, 1],
            z=particles[:, 2],
            px=particles[:, 3],
            py=particles[:, 4],
            pz=particles[:, 5],
            clock=particles[:, 6] * NS_TO_S,  # ns → s
            charge=particles[:, 7],
            status=particles[:, 8].astype(np.int32),
            ref_time_ns=float(header[0]),
            ref_momentum_eVc=float(header[1]),
            total_charge_nC=float(header[2]),
            ref_x_m=float(header[3]) if len(header) > 3 else 0.0,
            ref_y_m=float(header[4]) if len(header) > 4 else 0.0,
            source=str(path),
            format="astra_binary",
        )

        logger.info(
            f"Read ASTRA binary '{path.name}': "
            f"p_ref={header[1] * 1e-6:.2f} MeV/c, "
            f"Q={dist.total_charge_nC:.4f} nC, "
            f"N={dist.n_particle} ({dist.n_active} active)"
        )

        return dist

    def _read_ascii(self, path: Path) -> Distribution:
        """Read ASCII ASTRA format (space-separated values).

        Supports two variants:
          - With 5-value header line: ref_time ref_energy Q_total ref_x ref_y
          - Without header: particle data starts on line 1 (Generator format)

        Also normalizes Generator status conventions:
          Generator uses status=1 for active particles (unlike ASTRA's status=0).
        """
        with open(path) as f:
            lines = [line.strip() for line in f if line.strip() and not line.startswith("#")]

        if not lines:
            raise ValueError(f"Empty ASTRA ASCII file: '{path}'")

        # ── Detect header line ──
        # A header line has exactly 5 values. A particle line has 9-10 values.
        first_parts = lines[0].split()
        has_header = len(first_parts) == 5
        data_start = 1 if has_header else 0

        header = np.zeros(5)
        if has_header:
            for i, val in enumerate(first_parts[:5]):
                header[i] = float(val)
            logger.info(f"Detected 5-value header in '{path.name}'")

        # ── Parse particles ──
        data_lines = lines[data_start:]
        n_cols = len(data_lines[0].split())
        n_particles = len(data_lines)

        particles = np.zeros((n_particles, 9))
        for i, line in enumerate(data_lines):
            vals = line.split()
            for j in range(min(len(vals), 9)):
                particles[i, j] = float(vals[j])

        # ── Normalize Generator status conventions ──
        # Generator output: status=1 means active (in the bunch).
        # ASTRA tracking: status=0 means active.
        # If no status=0 particles are found, remap the most common status to 0.
        raw_status = particles[:, 8]
        if np.all(raw_status != 0):
            unique_vals, counts = np.unique(raw_status, return_counts=True)
            most_common = unique_vals[np.argmax(counts)]
            particles[raw_status == most_common, 8] = 0.0
            logger.info(
                f"Generator format detected: remapped status={int(most_common)} → 0 "
                f"(active) for {int(counts[np.argmax(counts)])} particles"
            )

        dist = Distribution(
            x=particles[:, 0],
            y=particles[:, 1],
            z=particles[:, 2],
            px=particles[:, 3],
            py=particles[:, 4],
            pz=particles[:, 5],
            clock=particles[:, 6] * NS_TO_S,
            charge=particles[:, 7],
            status=particles[:, 8].astype(np.int32),
            ref_time_ns=float(header[0]),
            ref_momentum_eVc=float(header[1]),
            total_charge_nC=float(header[2]),
            ref_x_m=float(header[3]) if len(header) > 3 else 0.0,
            ref_y_m=float(header[4]) if len(header) > 4 else 0.0,
            source=str(path),
            format="astra_ascii",
        )

        logger.info(f"Read ASTRA ASCII '{path.name}': N={dist.n_particle} ({dist.n_active} active)")

        return dist

    @staticmethod
    def _looks_like_index(col: np.ndarray) -> bool:
        """Check if a column looks like sequential particle indices."""
        if len(col) < 2:
            return False
        # Check if values are roughly sequential integers starting near 1
        diffs = np.diff(col)
        # Particle indices should be approximately 1.0 apart
        return bool(np.all(np.abs(diffs - 1.0) < 0.1))

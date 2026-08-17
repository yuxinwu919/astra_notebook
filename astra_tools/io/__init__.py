"""I/O layer: ASTRA distribution and evolution-file readers.

批 3: 公开面收敛 — 全部读取器在此统一 re-export, 前端不再
用 __import__ 深挖内部模块。
"""

from .astra_dist import AstraDistributionReader, write_distribution
from .astra_emit import (read_emit_files, read_ref_file, read_sigma_file,
                         read_cemit_file, read_log_file, parse_output_file,
                         output_file_type, EmitData, EmitSet, RefData,
                         SigmaData)
from .astra_misc import (read_pscan, read_scan, read_error, read_lab_file,
                         read_track_file, read_cathode_file, read_xemit2,
                         read_tremit, read_cr_emit, read_larmor, read_density,
                         read_tcheck)
from .field_map import (read_cavity_field, read_solenoid_field,
                        read_wake_potential, read_3d_field_map,
                        read_3d_field_map_components, FieldMap3D,
                        parse_field_map_file, expand_tws_field_map,
                        fix_laser_map_header, CavityField, SolenoidField,
                        WakePotential, TEField, read_te_field)
from .plot_steering import read_plot_steering, cp_index_colors


def read_distribution(path, fmt: str = "auto"):
    """Read a particle distribution file, auto-detecting the format.

    Currently only the ASTRA format (binary or ASCII) is supported.

    Args:
        path: path to the distribution file.
        fmt: format hint; 'auto' (default) or 'astra'.

    Returns:
        Distribution instance.
    """
    from pathlib import Path

    path = Path(path)
    if fmt not in ("auto", "astra"):
        raise ValueError("unknown format " + repr(fmt))
    reader = AstraDistributionReader()
    if fmt == "astra" or reader.probe(path):
        return reader.read(path)
    raise ValueError("cannot determine format of " + str(path))


__all__ = [
    "read_distribution", "AstraDistributionReader", "write_distribution",
    "read_emit_files", "read_ref_file", "read_sigma_file",
    "read_cemit_file", "read_log_file", "parse_output_file",
    "output_file_type", "EmitData", "EmitSet", "RefData", "SigmaData",
    "read_pscan", "read_scan", "read_error", "read_lab_file",
    "read_track_file", "read_cathode_file", "read_xemit2", "read_tremit",
    "read_cr_emit", "read_larmor", "read_density", "read_tcheck",
    "read_cavity_field", "read_solenoid_field", "read_wake_potential",
    "read_3d_field_map", "parse_field_map_file", "expand_tws_field_map",
    "read_3d_field_map_components", "FieldMap3D",
    "fix_laser_map_header", "CavityField", "SolenoidField",
    "WakePotential", "TEField", "read_te_field",
    "read_plot_steering", "cp_index_colors",
]

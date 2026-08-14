"""Namelist I/O for ASTRA & Generator input files.

write_namelist  - Python dict -> Fortran namelist block
format_namelist - normalize/format an existing input file
parse_namelist  - namelist block -> Python dict
"""

from .write import write_namelist
from .format import format_input_text, format_input_file
from .parse import parse_namelists

__all__ = ["write_namelist", "format_input_text", "format_input_file", "parse_namelists"]

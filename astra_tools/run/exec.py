"""Locate and run the ASTRA / Generator executables.

Search order for `name` (e.g. 'astra', 'generator'):
  1. system PATH
  2. <project_dir>/astra/<name>   (README layout)
  3. <project_dir>/ASTRA/<name>   (case-insensitive fallback; the repo
     actually uses the uppercase directory - fixed from the legacy
     code, which only probed lowercase 'astra/' and failed on
     case-sensitive filesystems)
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from pathlib import Path
from typing import Optional

logger = logging.getLogger(__name__)


def check_executable(name: str, project_dir: Optional[Path] = None) -> str:
    """Return the full path of an executable or raise FileNotFoundError."""
    found = shutil.which(name)
    if found:
        logger.info("found %s on PATH: %s", name, found)
        return found

    if project_dir is not None:
        project_dir = Path(project_dir)
        candidates = []
        for dirname in ("astra", "ASTRA"):
            p = project_dir / dirname / name
            if p.exists():
                candidates.append(p)
        for p in candidates:
            if os.access(p, os.X_OK) or p.suffix == ".exe":
                logger.info("found %s in project dir: %s", name, p)
                return str(p)
        # On Windows the extension may differ
        for p in candidates:
            if p.exists():
                logger.warning("%s exists but is not executable: %s", name, p)
                return str(p)

    raise FileNotFoundError(
        "executable '%s' not found. Put it on PATH or into "
        "<project>/astra/ or <project>/ASTRA/." % name
    )


def get_version(exe_path: str, work_dir: Path, timeout: int = 5) -> str:
    """Run the executable and extract a version string from its output."""
    try:
        result = subprocess.run(
            [exe_path], capture_output=True, text=True, timeout=timeout,
            cwd=str(work_dir), input="",
        )
    except (subprocess.TimeoutExpired, OSError) as e:
        logger.warning("version detection failed: %s", e)
        return "无法检测"
    for line in result.stdout.splitlines():
        if "version" in line.lower():
            return line.strip()
    return "未知版本"


def run_program(
    exe_path: str,
    work_dir: Path,
    input_file: Optional[str] = None,
    timeout: int = 3600,
    stream: bool = True,
) -> subprocess.CompletedProcess:
    """Run ASTRA / Generator in a working directory.

    Args:
        exe_path: executable path.
        work_dir: working directory (input files live here).
        input_file: optional input deck file name (ASTRA reads it from
            stdin when omitted; Generator takes it as argument).
        timeout: wall-clock limit [s] (default 1 h; ASTRA runs can be long).
        stream: stream output live instead of buffering.

    Raises:
        RuntimeError: non-zero exit code or timeout.
    """
    exe_name = Path(exe_path).name
    cmd = [exe_path] + ([input_file] if input_file else [])
    logger.info("running %s in %s", " ".join(cmd), work_dir)

    try:
        if stream:
            proc = subprocess.Popen(
                cmd, cwd=str(work_dir), text=True,
                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
            )
            lines = []
            assert proc.stdout is not None
            for line in proc.stdout:
                line = line.rstrip("\n")
                print(line)
                lines.append(line)
            returncode = proc.wait(timeout=timeout)
            result = subprocess.CompletedProcess(cmd, returncode, "\n".join(lines), "")
        else:
            result = subprocess.run(
                cmd, cwd=str(work_dir), capture_output=True, text=True,
                timeout=timeout,
            )
    except subprocess.TimeoutExpired:
        raise RuntimeError("%s timed out after %d s" % (exe_name, timeout))

    if result.returncode != 0:
        raise RuntimeError(
            "%s returned exit code %d\n%s"
            % (exe_name, result.returncode, result.stderr or result.stdout)
        )
    logger.info("%s finished successfully", exe_name)
    return result


def discover_outputs(work_dir: Path, stem: str, run: str = "001"):
    """Discover the output files of one run in a directory.

    Groups ASTRA output files by their base name:
      <stem>.Xemit.<run>, <stem>.Yemit.<run>, ..., <stem>.<NNNN>.<run>

    Returns a dict: {'emit': ..., 'sigma': ..., 'ref': ..., 'log': ...,
                     'phase': [paths], 'dist': [paths]}
    """
    work_dir = Path(work_dir)
    out = {"emit": None, "sigma": None, "ref": None, "log": None,
           "phase": [], "dist": [], "cemit": None}
    for f in sorted(work_dir.glob("*")):
        parts = f.name.split(".")
        if parts[0] != stem:
            continue
        if len(parts) >= 3 and parts[-1] == run:
            ext = parts[1]
            if ext == "Xemit":
                out["emit"] = f
            elif ext == "Sigma":
                out["sigma"] = f
            elif ext == "ref":
                out["ref"] = f
            elif ext == "Log":
                out["log"] = f
            elif ext == "Cemit":
                out["cemit"] = f
            elif ext.isdigit():
                out["phase"].append(f)
        elif f.suffix == ".ini":
            out["dist"].append(f)
    out["phase"].sort()
    return out


def backup_directory(src: Path, backup_root: Path):
    """Copy a simulation directory to <backup_root>/<YYYYmmdd_HHMMSS>."""
    from datetime import datetime
    import shutil as _shutil

    src = Path(src)
    if not src.exists():
        raise FileNotFoundError("simulation directory not found: " + str(src))
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    dst = Path(backup_root) / timestamp
    dst.parent.mkdir(parents=True, exist_ok=True)
    _shutil.copytree(src, dst)
    logger.info("backed up %s -> %s", src, dst)
    return dst

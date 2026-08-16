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
import re
import shutil
import subprocess
import threading
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
        input_file: optional input deck file name (passed as command
            line argument). stdin is always /dev/null so a batch binary
            can never block waiting for console input.
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
                # Jupyter 内核的 stdin 是无人写入的阻塞管道, 若批处理
                # 程序意外读 stdin 会永久挂起 ("跑几分钟不结束")。
                # 一律给 /dev/null: 读 stdin 立即得到 EOF。
                stdin=subprocess.DEVNULL,
            )
            lines = []

            def _drain() -> None:
                assert proc.stdout is not None
                # 批 5: 只保留尾部 (长跑数万行时内存有界); 失败回放用尾部
                for line in proc.stdout:
                    line = line.rstrip("\n")
                    print(line)
                    lines.append(line)
                    if len(lines) > 2000:
                        del lines[:1000]

            reader = threading.Thread(target=_drain, daemon=True)
            reader.start()
            try:
                returncode = proc.wait(timeout=timeout)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait(timeout=10)
                raise RuntimeError(
                    "%s timed out after %d s" % (exe_name, timeout))
            reader.join(timeout=10)
            result = subprocess.CompletedProcess(cmd, returncode, "\n".join(lines), "")
        else:
            result = subprocess.run(
                cmd, cwd=str(work_dir), capture_output=True, text=True,
                timeout=timeout, stdin=subprocess.DEVNULL,
            )
    except subprocess.TimeoutExpired:
        raise RuntimeError("%s timed out after %d s" % (exe_name, timeout))

    combined = (result.stderr or "") + "\n" + (result.stdout or "")
    if result.returncode != 0:
        raise RuntimeError(
            "%s returned exit code %d\n%s" % (exe_name, result.returncode, combined)
        )
    # ASTRA 解析失败时可能仍以 0 退出, 检查失败标记。
    # 注: "ERROR" 只认行首独立出现, 避免 Head 等正常文本里含
    # "ERROR" 字样时误报 (例如 'ERROR STUDY' 作为标题)。
    # 批 4: 补崩溃类标记 (ASTRA 段错误时可能仍以 0 退出)
    for marker in ("Program stops", "Error reading", "Segmentation fault",
                   "SIGSEGV", "core dumped", "Abort trap", "Trace/BPT trap"):
        if marker in combined:
            raise RuntimeError("%s 失败 (%s):\n%s" % (exe_name, marker, combined[-2000:]))
    if re.search(r"(?m)^\s*ERROR\b", combined):
        raise RuntimeError("%s 失败 (ERROR):\n%s" % (exe_name, combined[-2000:]))
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
    out = {"emit": None, "yemit": None, "zemit": None, "sigma": None,
           "ref": None, "log": None, "landf": None,
           "phase": [], "dist": [], "cemit": None}
    for f in sorted(work_dir.glob("*")):
        parts = f.name.split(".")
        if parts[0] != stem:
            continue
        if len(parts) >= 3 and parts[-1] == run:
            ext = parts[1]
            if ext == "Xemit":
                out["emit"] = f
            elif ext == "Yemit":
                out["yemit"] = f
            elif ext == "Zemit":
                out["zemit"] = f
            elif ext == "Sigma":
                out["sigma"] = f
            elif ext == "ref":
                out["ref"] = f
            elif ext == "Log":
                out["log"] = f
            elif ext == "LandF":
                out["landf"] = f
            elif ext == "Cemit":
                out["cemit"] = f
            elif ext.lstrip("-").isdigit():
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
    # 同秒内多次备份: 追加 -1/-2 后缀, 不覆盖不抛异常
    n = 0
    while dst.exists():
        n += 1
        dst = Path(backup_root) / (timestamp + "-" + str(n))
    dst.parent.mkdir(parents=True, exist_ok=True)
    _shutil.copytree(src, dst)
    logger.info("backed up %s -> %s", src, dst)
    return dst

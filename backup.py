#!/usr/bin/env python3
"""
备份 simulation_output 到 data/ 目录，以日期时间命名。

用法：
    python backup.py
"""

import shutil
from datetime import datetime
from pathlib import Path

PROJECT_DIR = Path(__file__).resolve().parent
SIM_DIR = "simulation_files"
BACKUP_ROOT = "data"


def main():
    src = PROJECT_DIR / SIM_DIR
    if not src.exists():
        print(f"✗ 仿真输出目录不存在: {src}")
        return

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    dst = PROJECT_DIR / BACKUP_ROOT / timestamp
    dst.parent.mkdir(parents=True, exist_ok=True)

    shutil.copytree(src, dst)
    print(f"✓ 仿真结果已备份到: {dst}")


if __name__ == "__main__":
    main()

"""astra-notebook 引导: 定位后端包 astra_tools 并应用全局样式.

用法 (任一 notebook 第一个代码单元格):
    %run _bootstrap.py

或从项目根目录启动 Jupyter 时直接 import。本文件按 __file__ 定位
项目根, 因此整个项目文件夹复制到任意路径都能工作。
"""

import sys
from pathlib import Path

# 项目根: 向上找到包含 astra_tools 的目录 (notebook 可位于任意
# 深度, 例如 notebooks/ 或 examples/<算例>/ 下)
_p = Path(__file__).resolve().parent
while not (_p / "astra_tools").is_dir() and _p != _p.parent:
    _p = _p.parent
PROJECT_ROOT = _p
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

import astra_tools
from astra_tools.plot.style import set_style

set_style()

# 默认工作目录 (运行时生成, 不进 git)
SIM_DIR = PROJECT_ROOT / "data" / "workspace"
SIM_DIR.mkdir(parents=True, exist_ok=True)

print("astra-notebook 后端已加载 (v%s)" % astra_tools.__version__)
print("项目根目录: %s" % PROJECT_ROOT)
print("模拟工作目录: %s" % SIM_DIR)

# 定位可执行文件 (PATH -> 项目 ASTRA/ 目录); 未找到时置 None 并告警
from astra_tools.run import check_executable

try:
    ASTRA_EXE = check_executable("astra", project_dir=PROJECT_ROOT)
    print("ASTRA    :", ASTRA_EXE)
except FileNotFoundError as _e:
    ASTRA_EXE = None
    print("警告: 未找到 astra 可执行文件 (PATH 或项目 ASTRA/ 目录)")
try:
    GENERATOR_EXE = check_executable("generator", project_dir=PROJECT_ROOT)
    print("Generator:", GENERATOR_EXE)
except FileNotFoundError as _e:
    GENERATOR_EXE = None
    print("警告: 未找到 generator 可执行文件")

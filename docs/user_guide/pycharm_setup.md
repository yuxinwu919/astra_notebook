# 在 PyCharm 中打开本项目

## 现状 (已配置好, 开箱即用)

* 虚拟环境 .venv/ 就在项目根目录 (Python 3.14.6, 全部依赖已装);
* requirements.txt 在根目录 — PyCharm 打开项目会自动识别:
  若提示 Install requirements from requirements.txt 可选择安装
  (用自带的 .venv 则跳过);
* Jupyter 内核 astra-notebook 已注册到 .venv
  (路径 .venv/share/jupyter/kernels/astra-notebook);
* ASTRA / Generator 可执行文件在 PATH
  (/Users/yuxinwu/programs/ASTRA/), 后端可自动定位;
* .idea/ 已在 .gitignore, IDE 配置不入库。

## 打开步骤

1. PyCharm → Open → 选择项目根目录 astra_notebook;
2. 解释器: Settings → Project → Python Interpreter →
   Add Interpreter → Existing → 选 .venv/bin/python
   (PyCharm 通常自动检测到 .venv 并列出);
3. (可选) 若弹出 requirements 安装提示且想用全新环境: 点 Install,
   完成后用 python -m ipykernel install --prefix .venv --name astra-notebook
   注册内核; 用自带 .venv 则跳过;
4. 打开 notebooks/ 下任意 .ipynb: 内核选 astra-notebook;
   按 01 → 02 → 03/04/05 顺序运行;
5. 运行测试: 终端执行
   .venv/bin/python -m pytest test/ -q  (357 项) 或
   在 PyCharm 的 pytest 配置里选 .venv 解释器、目录 test/。

## 目录速览

    notebooks/    5 个任务式 notebook (前端, 按编号使用)
    examples/     官方算例输入/输出 + 示例 notebook + 功能演示 demo
    astra_tools/  纯 Python 后端包 (复制即用, 无打包)
    test/         五层测试 (见 docs/dev_manual/test_plan.md)
    docs/         用户手册 / 开发手册 / 物理备忘录
    data/         本地运行产物 (gitignored)

## 常见问题

* notebook 报 kernel not found: Settings → Jupyter → 勾选
  Use Jupyter server from project interpreter, 内核列表里选
  astra-notebook; 仍无则重跑 ipykernel 注册命令。
* astra / generator 找不到: 检查 Settings → Tools → Terminal
  是否继承系统 PATH (默认继承); 或把可执行文件复制到项目
  ASTRA/ 目录 (后端第二搜索位置)。
* matplotlib 字体缓存警告: 首次绘图会构建字体缓存, 属正常。
* 若系统同时装了 anaconda 且 shell 里激活过 conda: Jupyter 内核列表
  可能混入 anaconda 的 python3, 且 python -m jupyter 会误用其
  site-packages。对策: 内核手动选 astra-notebook; 跑测试用仓库的
  bash test/e2e_notebooks.sh (脚本已内置 JUPYTER_PATH 隔离, 不受
  conda 影响)。

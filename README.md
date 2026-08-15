# astra-notebook

**ASTRA / Generator 前后处理工作台** — 以 Jupyter Notebook 为前端、Python 包为后端,
替代 macOS 上缺失的官方 postpro / lineplot / fieldplot 图形程序。

- **前端**: 8 个任务式 Notebook, 参数用表单点选 (元数据驱动, 覆盖 ASTRA 手册第 6 章
  全部 13 个 namelist 与第 7 章 Generator 全部参数), 也可直接使用现成 .in 文本;
  结果以统计表格 + 现代 KDE 密度图呈现, 支持 CSV/npz 数据导出。
- **后端**: astra_tools 纯 Python 包 (不打包、不分发), **复制整个文件夹即可使用**;
  统计与绘图代码全部经物理审查, 与 ASTRA 自身输出交叉验证 (误差 < 0.02%)。
- **覆盖**: postpro/lineplot/fieldplot 三大图形程序的功能替代 — 相空间与切片、
  全部演化曲线 (z 或 t 轴)、发射度/能散/光学函数、场图 (1D/TWS/3D 截面/轴上
  剖面)、孔径叠加、阴极发射、激光与等离子体、核心电荷分数曲线、
  BFF (直接法 + FFT 快速路径)、PScan/Scan/Error 扫描图。
- **验证**: 83 项测试 (五层, 见 docs/dev_manual/test_plan.md), 含真实 ASTRA
  PScan/Scan 交叉验证; 官方 9 算例一键复现 + 黄金比对 (08 号 notebook)。

## 快速开始

    # 1. 依赖 (Python 3.9+)
    python3 -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt
    python -m ipykernel install --prefix .venv --name astra-notebook

    # 2. ASTRA/Generator 可执行文件: 放入 PATH 或项目 ASTRA/ 目录

    # 3. 启动 Jupyter
    jupyter notebook

    # 4. 按编号运行 notebooks/ 下的 Notebook: 01 -> 02 -> 03/04/05

## Notebook 一览

| # | Notebook | 对应原程序 | 作用 |
|---|----------|-----------|------|
| 01 | generator | generator | 初始束团生成 (参数表单 -> 运行 -> 束团预览) |
| 02 | astra | astra | 工作区自检 + 追踪设置 (全部 namelist 表单 + .in 文本) + 运行/输出清单 |
| 03 | postpro | postpro | 相空间分析/统计/切片/BFF/导出 (含孔径叠加与核心曲线) |
| 04 | lineplot | lineplot | 束流参数演化 (九图 + 速度/步长 + t 轴变体) |
| 05 | fieldplot | fieldplot | 腔场/螺线管/3D 场图/激光/等离子体 |
| 06 | examples | — | 教程与功能总览: 六章按 01→05 逐章教学并展示全部功能 (真实算例 + 合成教学数据) + 一键复现与黄金比对 |

## 目录结构

    astra_notebook/           # 复制这一层即可使用
    ├── astra_tools/          # 后端包 (io/analysis/plot/deck/run/widgets/export)
    ├── notebooks/            # 6 个前端 Notebook + _bootstrap.py
    ├── examples/             # DESY 官方算例 (输入 + 黄金输出)
    ├── test/                 # 五层测试
    ├── docs/                 # user_guide / dev_manual / physics_notes
    ├── data/                 # 运行时工作区 (gitignore)
    └── AGENTS.md             # 协作规范 (英文, 给后续开发 agent)

## 物理约定 (全部经真实 ASTRA 输出验证, 详见 docs/physics_notes/)

- 发射度单位 **π mm mrad**: 数值上即 eps_rms 的 mm·mrad (π 表示 RMS 相空间椭圆面积
  语义), 与 ASTRA 打印值完全一致; 换算 SI = ×1e-6, 不乘 π。
- 活跃粒子: status > 1 (手册 4.13); z/pz/clock 在文件中相对参考粒子, 读取时
  归一化为绝对坐标; 螺线管场中发射度用正则动量 p̃x = px + c·Bz·y/2。

## 测试

    .venv/bin/python -m pytest test/ -q     # 单元 + 黄金回归 + 交叉验证 + 绘图正确性
    # 端到端 (需可执行文件): 逐个执行 notebooks/ 下的 8 个 Notebook

## 文档

- 用户手册: docs/user_guide/ (中文; 含 PyCharm 打开指南 pycharm_setup.md)
- 开发手册: docs/dev_manual/ (中文; 含完整测试方案 test_plan.md 与
  可委托 Codex 的任务清单 codex_tasks.md)
- 物理审查备忘录: docs/physics_notes/ (约定 01-09 + 未解决物理问题报告 10)
- 环境: 自带 .venv (Python 3.14.6, 依赖已装); requirements.txt 供
  PyCharm 自动识别一键安装

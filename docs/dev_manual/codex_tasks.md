# 任务分工: 主代理自留 vs Codex (PyCharm MCP)

> 2026-08-15 更新。原则: **我做得了、做得好的一律自己做;
> 只有真正需要 IDE/PyCharm 的任务才交给 Codex**。
> Codex 通过 MCP 连接 PyCharm, 但后端的一切物理/工程验证
> (跑 ASTRA、算统计、写测试) 主代理用 bash+python 完全可做,
> 不需要 IDE。

## 环境接入 (Codex 使用前)

当前状态 (2026-08): 115 项测试、e2e 18 本 (5 任务式 + 13 examples/)。
已完成: 示例 notebook 拆分与平铺、postpro 步进器 + 演示 notebook、
手册 5.5-5.7 覆盖审计 (coverage_audit.md)、散点渲染替换 KDE、
4 本功能演示 demo (generator_demo / bff_demo /
stats_validation_demo / lineplot_demo)。

1. PyCharm 打开项目根目录 astra_notebook;
2. 解释器选 .venv/bin/python (Python 3.14.6);
3. 测试: .venv/bin/python -m pytest test/ -q;
4. 依赖: requirements.txt (PyCharm 会提示安装);
5. 铁律见 AGENTS.md (物理 8 条 + 测试五层 + 不做打包);
   改物理约定必须先改测试并给手册依据。

---

# 第一部分: 主代理自留任务 (物理/工程验证, 按优先级)

## M1: read_error 合成单元测试 — **已完成** (test/test_misc_readers.py)

## M2: lab.001 标签读取器 — **已完成** (read_lab_file + Scan/Error
图标题; golden: examples/Manual_Example/golden/Example.lab.001)

## M3: Cemit 核心发射度真实交叉验证 (高, 3-4h)
physics_notes/10 问题 2。Manual_Example 真跑 C_EmitS=T ->
golden Cemit.001; 用 analysis/core.py 的 0.8/0.9/0.95 核心分数曲线
逐列对照 C80/90/95; 不吻合则试横向核心/相空间距离核心 (手册
4.13.8); 结论写 physics_notes/10。

## M4: Sigma 3.83 因子 — **已解决** (2026-08)
文件把动量列归一化到 mc、能量列归一化到 mc^2; 3.83 = 1/(mc[MeV])^2。
read_sigma_file 已换算 SI, 归一化 eigen-emittance 与 Xemit 对照 <8%,
enz 与 Zemit <0.01% (physics_notes/06, test_cross_validation.py)。

## M5: Error 真实运行 + golden (中, 2-3h)
physics_notes/10 问题 4。Manual_Example 加 &ERROR (ErrorS=T,
Err_MaxB) 真跑 -> Error.001 golden; 与 Scan 名义值对照。

## M6: Plasma_2 激光尾场物理量级验证 (中, 2-3h)
physics_notes/10 问题 8。按手册第 8 章尾场公式估算 a0~0.8、
1e17 cm^-3、800nm 的能增量级, 对照 golden 末态 421.8 MeV;
检查 laser.dat 单位约定 (a0² vs a0)。

## M7: slice 发射度临时第三方对照 (一次性, 不进仓库)
physics_notes/10 问题 3。用 venv 里的 pmd-beamphysics 对
Example.0150.001 做 slice_analysis, 与 compute_slice_analysis 对照;
结论写入 physics_notes/10。(铁律: 不进仓库、不成为运行时依赖)

## M8: 持续回归
每次改动: 红测试先行 -> pyflakes 清 -> 全量 pytest ->
涉及 notebook 重跑 e2e (bash test/e2e_notebooks.sh)。

---

# 第二部分: Codex (PyCharm MCP) 专属任务 (只有 IDE 做得好)

## C1: 5 本任务式 + examples/ 教学/演示 notebook 的交互 UX 审阅 (推荐先做)

* 背景: notebook 刚精简为 01_generator / 02_astra / 03_postpro /
  04_lineplot / 05_fieldplot 以及 examples/ 下 13 本教学/演示
  notebook (含 postpro_step_demo.ipynb 与 4 本功能演示 demo),
  从未做过人工交互审阅。
* 步骤: 在 PyCharm 的 Jupyter 面板逐个打开运行 (内核
  astra-notebook), 检查: widget 表单渲染与默认值 (02_astra)、
  下拉选择器选项 (03_postpro)、图表尺寸/图例/中文文案是否清晰、
  暗色主题下的可读性、报错信息是否指向正确 notebook;
  对每个问题截图并给出定位 (cell 编号)。
* 输出: 问题清单 (截图 + cell 定位 + 建议), 提交到主代理修复。
* 注意: 只读审阅, 不改代码。

## C2: PyCharm 工程配置 (一次性)

* 创建并验证运行配置: (a) pytest: 解释器 .venv, 目录 test/;
  (b) e2e: Shell Script test/e2e_notebooks.sh;
  (c) Jupyter Server: 内核 astra-notebook, 工作目录项目根;
* 配置 ruff (或保留 pyflakes) 作为外部工具, 命令:
  .venv/bin/python -m pyflakes astra_tools test;
* 验证 02_astra 的 namelist_form 在 PyCharm Jupyter 里可交互修改;
* 输出: .idea/ 里的配置说明 (不入库, 写 docs/user_guide/pycharm_setup.md 增补)。

## C3: IDE 静态检查告警清单 (可选, 定期)

* 运行 PyCharm Inspections (全项目), 过滤风格类告警,
  只输出: 未解析引用 / 可能为 None 的解引用 / 类型不匹配 /
  未使用参数; 提交清单供主代理逐一核实修复。

## C4: laser.dat 读取性能/内存剖析 (可选)

* 用 PyCharm profiler 对 read_3d_field_map('examples/Plasma_Example_2/laser.dat',
  81×81×400, 65MB) 剖析: 当前实现逐 token float() 解析, 记录耗时与峰值内存,
  若 > 30s 建议 np.fromstring 加速方案 (只提方案, 不改代码)。

---

## 通用提交要求 (双方)

* 红测试先行 -> 实现 -> pyflakes 清 -> 全量 pytest -> 涉及 notebook
  重跑 e2e -> 单条中文提交说明 (含验证了什么)。
* 不改 AGENTS.md 物理 8 条、不改 golden 数据 (除非任务明确要求)。

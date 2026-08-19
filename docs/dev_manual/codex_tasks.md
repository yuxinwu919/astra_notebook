# 任务分工: 主代理自留 vs Codex (PyCharm MCP)

> 2026-08-15 更新。原则: **我做得了、做得好的一律自己做;
> 只有真正需要 IDE/PyCharm 的任务才交给 Codex**。
> Codex 通过 MCP 连接 PyCharm, 但后端的一切物理/工程验证
> (跑 ASTRA、算统计、写测试) 主代理用 bash+python 完全可做,
> 不需要 IDE。

## 环境接入 (Codex 使用前)

当前状态 (2026-08): 357 项测试、e2e 19 本 (5 任务式 + 14 examples/)。
已完成: 示例 notebook 拆分与平铺、postpro 步进器 + 演示 notebook、
手册 5.5-5.7 覆盖审计 (coverage_audit.md)、散点渲染替换 KDE、
6 本功能演示 demo (generator_demo / bff_demo /
stats_validation_demo / lineplot_demo / fieldplot_demo / postpro_demo)。

1. PyCharm 打开项目根目录 astra_notebook;
2. 解释器选 .venv/bin/python (Python 3.14.6);
3. 测试: .venv/bin/python -m pytest test/ -q;
4. 依赖: requirements.txt (PyCharm 会提示安装);
5. 铁律见 AGENTS.md (物理 8 条 + 测试五层 + 不做打包);
   改物理约定必须先改测试并给手册依据。

---

# 第三阶段路线 (2026-08-19, 两轮对抗性审计+修复后)

> 两轮审计已修复 18 个 P1 与 ~30 个 P2/P3, 测试 357→441 全绿
> (commit 87e7e96)。本节规划补齐审计暴露的**验证空缺**, 按优先级:

## R1 (P0 级验证缺口) 真实二进制分布 golden — 2h
ASTRA 默认输出即 binary, 但仓库仍无真实二进制文件对照。Manual_Example
设 BINARY=T 真跑 -> 归档 golden (首个二进制分布文件); 与
AstraDistributionReader 逐字节对照 (消歧/字节序/相对坐标/头 Q=|Q|),
并把 golden 挂进 test_cross_validation。验证 2026-08 重写的 9/10 列
语义判定与种类列解析在真实文件上成立。

## R2 (P0 级验证缺口) 低能/大发散交叉验证 — 3h
现有 golden 全部 ~1 GeV 低发散均匀电荷; 动能全动量 (|p|)、γβ 口径
等因子类修复在低能区未与 ASTRA 对照。做两个真跑:
  (a) 低能束 (generator 发射段或 ~MeV deck) -> Xemit/Zemit 对照;
  (b) 90deg 大发散探针束 -> Zemit 对照 (直接验证 pz-only→|p| 修复
      与 ASTRA 的口径一致; 若 ASTRA 自身用 pz-only 需在笔记 07 记录)。
结论写 physics_notes/07 与 10。

## R3 正电子/离子全链验证 — 2h
构造 species=2 分布 (index 列), 验证 canonical_signs 在统计/slice/
核心分数/任意相空间/相空间图 5 条路径的口径一致 (仓库无正电子算例);
纯内部解析验证即可 (无 ASTRA 正电子对照需求), 加回归测试。

## R4 Cemit 核心发射度真实对照 — 3-4h (原 M3, 高)
见 M3。README 已限定 +5~25% 未对齐; 真跑 C_EmitS=T 逐列对照,
不吻合则试横向核心/相空间距离核心定义 (手册 4.13.5)。

## R5 Error 文件真实 golden — 2-3h (原 M5)
ErrorS=T 真跑 -> Error.001 golden, 验证 read_error 列语义。

## R6 slice 第三方一次性对照 — 1h (原 M7, 不进仓库)
pmd-beamphysics 临时脚本对照 compute_slice_analysis。

## R7 e2e 全套终跑 — 2h
19 本 notebook (5 任务式 + 14 examples) 用本地 ASTRA 二进制执行
test/e2e_notebooks.sh, 确认两轮改动后全链 (含 postpro 的
_bunch_center_z/_bz_t 传播、lineplot 的 run_key) 无回归。

## R8 收尾卫生 — 1h
pyflakes 历史遗留 (~10 条未用导入/变量) 清理; C1 (Codex/PyCharm
交互 UX 审阅) 视人力安排; 完成后更新 physics_notes/README 状态表。

**执行纪律** (沿用 M8): 红测试先行 -> pyflakes 清 -> 全量 pytest ->
涉及 notebook 重跑 e2e; 改物理约定必须先改测试并给手册依据。
**预计总量**: ~16-18h 主代理工作量, 无外部阻塞 (ASTRA 二进制本地已有)。

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
文件把动量列与能量列都归一化到 mc (交叉项乘 mc 一次), 即
scale=[1, mc, 1, mc, 1, mc]; 3.83 = 1/(mc[MeV])^2 仅是 MeV 单位表象。
read_sigma_file 已换算 SI, 归一化 eigen-emittance 与 Xemit 对照 <8%,
enz 与 Zemit <0.01% (physics_notes/06, test_cross_validation.py)。

## M5: Error 真实运行 + golden (中, 2-3h)
physics_notes/10 问题 4。Manual_Example 加 &ERROR (ErrorS=T,
Err_MaxB) 真跑 -> Error.001 golden; 与 Scan 名义值对照。

## M6: Plasma_2 激光尾场物理量级验证 — **已完成** (2026-08, T9)
physics_notes/10 问题 8 已关闭。按 Astra-Manual V3.2 §6.9 线性激光
尾场公式独立积分 (a0=0.8, n=1e17 cm^-3, λ=800 nm, σz=6.4 μm,
w0=50 μm): 末态动能 422.7 MeV vs golden 421.8 MeV, 吻合 0.2%;
E_peak = 4.32 GV/m = 0.142·E_WB (线性区); 失相长度 1.84 m 与
泵浦耗尽 ~1.7 m 均 >> 0.1 m 靶长; 激光 ~54 TW << 自聚焦临界
~0.3 PW; 束加载可忽略; √(2π) 系数变体 (743.6 MeV) 被排除。
限制: DESY 官方说明 PDF 无期望末态能量数字, 以手册公式为准。

## M7: slice 发射度临时第三方对照 (一次性, 不进仓库)
physics_notes/10 问题 3。用 venv 里的 pmd-beamphysics 对
Example.0150.001 做 slice_analysis, 与 compute_slice_analysis 对照;
结论写入 physics_notes/10。(铁律: 不进仓库、不成为运行时依赖)

## M8: 持续回归
每次改动: 红测试先行 -> pyflakes 清 -> 全量 pytest ->
涉及 notebook 重跑 e2e (bash test/e2e_notebooks.sh)。

---

# 第二部分: Codex (PyCharm MCP) 专属任务 (只有 IDE 做得好)

## C1: 5 本任务式 + examples/ 示例/演示 notebook 的交互 UX 审阅 (推荐先做)

* 背景: notebook 为 generator / astra / postpro / lineplot /
  fieldplot 以及 examples/ 下 14 本示例/演示
  notebook (含 6 本功能演示 demo),
  从未做过人工交互审阅。
* 步骤: 在 PyCharm 的 Jupyter 面板逐个打开运行 (内核
  astra-notebook), 检查: widget 表单渲染与默认值 (astra)、
  下拉选择器选项 (postpro)、图表尺寸/图例/中文文案是否清晰、
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
* 验证 astra 的 namelist_form 在 PyCharm Jupyter 里可交互修改;
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

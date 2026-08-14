# Codex (IDE MCP) 任务清单 — 当前状态版

> 2026-08-15, 批 A/B/C 已由主代理完成 (83 测试通过)。以下是按当前
> 代码状态整理的、可整块委托给 Codex 的具体任务。每条含: 背景 /
> 步骤 / 验收 / 涉及文件 / 依赖。

## 环境接入 (每个任务开始前)

1. 用 PyCharm 打开项目根目录 astra_notebook;
2. 解释器选 .venv/bin/python (Python 3.14.6, 依赖已装好);
3. 测试运行: .venv/bin/python -m pytest test/ -q
4. 改动后跑 pyflakes: .venv/bin/python -m pyflakes astra_tools test
5. ASTRA/Generator 在 /Users/yuxinwu/programs/ASTRA/ (PATH 已含);
   真跑在 data/ 下 (gitignored), golden 产物才进 examples/。
6. 铁律: 见 AGENTS.md (物理 8 条 + 测试五层 + 不做打包)。改物理
   约定前必须先改对应测试并说明手册依据。

---

## T1: read_error 合成单元测试 (低风险, 约 1 小时)

* 背景: astra_tools/io/astra_misc.py 的 read_error (run#, z,
  FOM(1..10), 12 列) 无任何测试; 真实 ErrorS 运行见 T5。
* 步骤: 合成一个 12 列 Error.001 (numpy.savetxt), 断言列语义
  (run 为 int, FOM 形状 (N,10), z 原样); 空文件/缺列报 ValueError。
* 验收: 新测试放入 test/, pytest 绿, pyflakes 清。
* 文件: astra_tools/io/astra_misc.py, test/。

## T2: lab.001 标签读取器 + 接入 Scan/Error 图标题 (低, 约 2 小时)

* 背景: &SCAN 运行时生成 <stem>.lab.001 (A80 三行: X轴/Y轴/标题);
  目前 plot_scan_fom / plot_error_hist 用默认标题。示例文件:
  data/scan_validate/Example.lab.001。
* 步骤: 新增 read_lab_file(path) -> dict(xlabel, ylabel, title);
  plot_scan_fom / plot_error_hist 增加 lab 可选参数, 有则用其标题/
  轴名; 导出到 astra_tools/plot/__init__; 测试用真实
  Example.lab.001 (先提交为 golden)。
* 验收: 单元测试 + 绘图标签断言。
* 文件: astra_misc.py, advanced_plots.py, plot/__init__.py, test/。

## T3: Cemit 核心发射度真实交叉验证 (高价值, 约 3-4 小时)

* 背景: 见 physics_notes/10 问题 2。Cemit C80/90/95 算法未复现;
  我们已有 analysis/core.py (纵向 |q| 分数核心)。
* 步骤:
  1. 复制 examples/Manual_Example 输入到 data/cemit_validate/,
     在 &OUTPUT 加 C_EmitS=T, 真跑 astra -> Example.Cemit.001;
  2. 用 read_distribution 读 Example.0150.001,
     compute_core_fraction_curves(fractions=(0.8,0.9,0.95)),
     与 Cemit 末行 Cx80/Cx90/Cx95 对照;
  3. 若不吻合, 依次试: 横向 x 分数核心 / 相空间距离核心 / 亚束团
     定义 (手册 4.13.8) 并记录哪个定义能复现;
  4. 结论写入 physics_notes/10 (更新问题 2), 黄金 Cemit.001 提交
     examples/Manual_Example/, 测试进 test/。
* 验收: 若复现成功 -> 对照测试 rel<5%; 若三种定义都不吻合 -> 文档
  记录失败矩阵, 保留只读展示状态并说明。
* 文件: analysis/core.py (可能), io/astra_emit.py, physics_notes/10, test/。

## T4: Sigma 3.83 因子实验 (高价值, 约 3 小时)

* 背景: 见 physics_notes/10 问题 1 与 physics_notes/06。
  sig(2,2)/sig(6,6) 与 <p~x^2>/sigmaE^2 差 ~3.83 因子, 手册未说明。
* 步骤:
  1. 读 examples/Manual_Example/Example.Sigma.001 全部行;
  2. 对每行计算 R22 = sig22 / (p_ref^2 * sigx'^2) 与
     R66 = sig66 / sigE^2 (sigx'/sigE 取同 z 行 Xemit/Zemit),
     得到 R(z);
  3. 检验候选因子: (bg)^2, g^2, g(1+b^2), g^2(1+b^2), 2g^2-1,
     (g^2+1)/2 ... 与 R(z) 随 gamma 变化曲线 (gamma 从 ~2 到
     ~1000) 做最小二乘, 锁定唯一假设; 必要时在 data/ 下用不同
     能量真跑 2-3 个 SigmaS 样本加大 gamma 跨度;
  4. 更新 read_sigma_file 注释与 physics_notes/06、10。
* 验收: 找到确切因子并把 eigen-emittance 去掉实验性标记 (或证明
  无一致因子, 文档化证据)。
* 文件: io/astra_emit.py, physics_notes/06+10, test/ (若修正)。

## T5: Error 真实运行 + golden (中, 约 2-3 小时)

* 背景: physics_notes/10 问题 4。
* 步骤: Manual_Example 副本加 &ERROR (ErrorS=T, Err_MaxB(1)=0.02,
  参照手册 4.12) 真跑 -> Example.Error.001 golden; 对照名义值:
  Error 运行 FOM 与 Scan 名义点 FOM 一致性 (误差扰动下均值应回
  名义值); 测试照 T1 扩展。
* 验收: golden + 对照测试; 若 &ERROR 语法有出入, 以手册 4.12 为准
  并记录。
* 文件: examples/Manual_Example/golden/, test/。

## T6: README / 文档特性更新 (低, 约 1 小时)

* 背景: 批 C 新增了 t 轴变体、孔径叠加、激光/等离子体图、核心
  分数曲线、PScan/Scan 交叉验证。
* 步骤: 更新 README.md 的 Notebook 一览与功能列表; docs/user_guide/
  补相应说明; 保证 docs/dev_manual/README.md 与
  physics_notes/README.md 的索引包含 09/10 与 test_plan/codex_tasks。
* 验收: markdown 链接可点、无死链。

## T7: 待机/可选 — Plasma_Example_2 激光图 fixture (低)

* 背景: laser.dat 65MB 本地 gitignored; plot_laser_on_axis 已用
  3D_test.ex 验证。可选: 加一个 skipif(laser.dat 不存在) 的激光
  剖面测试, 供有该文件的机器使用。
* 文件: test/。

## T8: 待机/可选 — slice 发射度临时第三方对照 (中, 一次性)

* 背景: physics_notes/10 问题 3。ASTRA 无 slice 级输出。
* 步骤 (仅本地, 不进仓库): 临时脚本用 venv 里已有的
  pmd-beamphysics 的 slice_analysis 对 Example.0150.001 求 slice
  epsn, 与 compute_slice_analysis 对照; 结论写入 physics_notes/10。
* 注意: 铁律 pmd-beamphysics 不得成为运行时依赖 — 只允许测试期
  临时对照脚本, 不得被仓库代码 import。

## T9: Plasma_Example_2 激光尾场物理验证 (中, 约 2-3 小时)

* 背景: physics_notes/10 问题 8。laser.dat 图头在新构建下只能用
  整数计数形式; 已据此重新生成 golden, 但尾场强度的绝对正确性
  未独立验证 (末态 E 100->421.8 MeV)。
* 步骤: 按手册第 8 章 (场图与等离子体) 的尾场公式估算 1e17 cm^-3、
  a0~0.8、800 nm 驱动激光在 zeta=-45 um 处的轴向场与能增, 与
  421.8 MeV 对照 (量级一致即通过); 或查 DESY 示例页的参考结果;
  结论写入 physics_notes/10 问题 8。
* 验收: 量级论证写入备忘录; 若发现量级不符, 检查 laser.dat 的
  单位约定 (a0^2 vs a0) 并给出修正建议。
* 文件: physics_notes/10。

## 通用提交要求

* 每条任务: 红测试先行 -> 实现 -> pyflakes 清 -> 完整 pytest ->
  涉及 notebook 则重跑 e2e -> 单条中文提交说明 (含验证了什么)。
* 不改 AGENTS.md 物理 8 条、不改 golden 数据 (除非任务明确要求)。

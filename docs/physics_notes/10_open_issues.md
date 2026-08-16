# 10 — 未解决的物理问题报告

> 状态: 2026-08 更新。118 项测试、e2e 18 本通过; Sigma(问题1)已解决。
> 每条给出: 现象 / 已知证据 / 影响范围 / 建议实验 / 风险评级。

## 1. Sigma 矩阵 ~3.83 因子 — **已解决** (2026-08)

**结论**: 文件把动量列归一化到 mc、能量列归一化到 mc^2 (无量纲);
3.83 = 1/(mc[MeV])^2。read_sigma_file 现已换算到 SI 并导出归一化
eigen-emittance, 与 Xemit eps_n 对照 < 8% (enz 与 Zemit eps_zn
< 0.01%)。详见 physics_notes/06。

## 2. Cemit 95/90/80 核心发射度算法未独立复现 ★高

**现象**: Cemit 文件的 Cx95/Cx90/Cx80 列 (核心发射度) 目前只读、
只显示, 未用粒子数据独立复现其"取核心"的算法。

**已知证据**: 我们新写的核心机制 (analysis/core.py, 按纵向 |q| 分数
取中心) 输出与全束团统计在 fraction=1.0 严格一致, 但与 ASTRA 的
C80/90/95 尚无对照。

**建议实验**: Manual_Example 加 C_EmitS=T 真跑 → golden Cemit.001;
用 compute_core_fraction_curves(f=0.8/0.9/0.95) 逐列对照, 若不吻合
再试横向核心、相空间距离核心等定义。细节见 codex_tasks.md T3。
**风险**: 中 — Cemit 仅展示用。

## 3. slice 发射度与失配参数缺乏 ASTRA 级对照 ★中

**现象**: slice 发射度 (analysis/slices.py) 只有解析验证 (高斯束团
slice εn = σxσpx/mc²) 与 ζ≥1 下界 (失配参数); ASTRA 自身不输出
slice 级发射度文件, 无法做逐 slice 交叉验证。

**已知证据**: 批 A 已统一约定 (p_ref + 正则动量 + ddof=0); 解析
高斯测试 rel<2%; ζ 图仅有数值下界检查。
**建议实验**: 若需更强验证, 可在测试环境临时用 pmd-beamphysics 的
slice_analysis 对同一相空间文件比对 (仅本地临时脚本, 不进入仓库,
不成为运行时依赖)。见 codex_tasks.md T8。
**风险**: 低-中 — slice 分析用于展示与失配诊断。

## 4. Error 文件 (ErrorS) 无真实 golden ★中

**现象**: read_error 读者已就绪但从未用真实 ErrorS 运行验证; 手册
Table 4 列语义 (run#, z, FOM(1..10)) 仅靠格式推断。

**建议实验**: Manual_Example 加 &ERROR Err_MaxB(1)=0.02 等 → 
Error.001 golden, 与名义值运行 (Scan 或 Xemit) 对照。见
codex_tasks.md T5。
**风险**: 低 — Error 功能目前无用户入口。

## 5. PScan 与 ref 粒子的 0.38% 偏差 (已解释, 非缺陷) ★低

**现象**: PScan 相角 0 的能量 vs 同一运行 ref.001 末态参考粒子动能
差 0.38%。

**结论**: PScan 是单粒子扫描且不含空间电荷; 分布参考粒子带 z/clock
偏移, 注入相位差约 2°。同一运行 Zemit 束团平均能量与 PScan 峰差
2.3% (空间电荷 + 相位偏移), 属预期物理。测试容差取 0.5%,
记录于 docs/physics_notes/09。

## 6. BFF 无 ASTRA 级对照 ★低

**现象**: BFF (|F(k)|²) 只有解析高斯对照 (exp(−(kσz)²)) 与 FFT/直接
法互证; ASTRA 不输出 BFF 文件, 无法直接对照。FFT 路径在深零点处
有绝对误差地板 (binning 相位误差 ≤ π/512), 已在测试容差中记录。

## 7. 依赖本地大文件的展示 ★低

**现象**: Plasma_Example_2 的 laser.dat (65MB) 与本地的 laser 轴上
剖面图依赖该文件; 未纳入 CI。plot_laser_on_axis 已用 3D_test.ex
验证通用性。

## 8. Plasma_Example_2 激光 3D 图头在不同 ASTRA 构建间不可移植 ★中

**现象** (2026-08 实测, macOS Apple Silicon 构建): DESY 原版
laser.dat 头三行是 (n, min, spacing) 且计数写成浮点 (8.1e+01):
* 原版浮点计数   -> ASTRA "Error while reading file: laser.dat"
* 逐值网格行 (n 后跟 n 个值) -> 读取 3D 图时 SIGSEGV
* 整数计数 + (n,min,spacing) -> 正常 (fix_laser_map_header 现行
  输出, 已改为此形式)

**处置**: Plasma_Example_2 的 golden (Xemit/Yemit/Zemit/Log/0011.001)
与 golden_expected.json 已用整数计数头的真实运行重新生成
(末态 E_kin 100 -> 421.8 MeV, eps_nx = 0.156 μm); 旧 golden
(239.6 MeV, eps_nx = 108 μm) 来自旧头形式的运行, 无法在当前构建
复现, 已作废。新的绝对物理 (与 DESY 算例说明的期望尾场强度是否
一致) 尚未独立验证 — 建议委托 Codex 对照手册第 8 章尾场公式或
DESY 示例页 (见 codex_tasks.md T9)。

## 建议优先级

T3 (Cemit) → T4 (Sigma 3.83) → T5 (Error) → T9 (Plasma_2 激光物理)
→ T8 (slice 临时对照)。详见 docs/dev_manual/codex_tasks.md。

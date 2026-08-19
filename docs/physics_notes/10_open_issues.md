# 10 — 未解决的物理问题报告

> 状态: 2026-08-19 第三阶段完成后更新。507 项测试、e2e 19 本全绿。
> 问题 1/2/4/8 已关闭; 问题 3 的第三方对照 (R6) 因 venv 无
> pmd-beamphysics 延后; 问题 6 (BFF) 维持"无 ASTRA 级对照"限定
> (ASTRA 不输出 BFF)。两轮对抗性审计的全部修复与第三阶段验证结论
> 见 docs/dev_manual/codex_tasks.md 第三阶段路线。
> 每条给出: 现象 / 已知证据 / 影响范围 / 建议实验 / 风险评级。

## 1. Sigma 矩阵 ~3.83 因子 — **已解决** (2026-08)

**结论**: 文件把动量列与能量列都归一化到 mc (无量纲; 交叉项乘 mc
一次), 即 scale=[1, mc, 1, mc, 1, mc]; 3.83 = 1/(mc[MeV])^2 仅是
MeV 单位表象 (能量列不除以 mc^2)。read_sigma_file 现已换算到 SI 并导出归一化
eigen-emittance, 与 Xemit eps_n 对照 < 8% (enz 与 Zemit eps_zn
< 0.01%)。详见 physics_notes/06。

## 2. Cemit 95/90/80 核心发射度算法未独立复现 — **已关闭** (2026-08-19, R4 真跑破译)

**结论**: C_EmitS=T 真跑生成 golden Cemit.001, 逐列对照后破译 ASTRA
算法: 核心发射度 = 按单粒子振幅 J_i 升序取前 f·N 个粒子的**平均振幅**
ΣJ/N (**不是**核心子集重算 rms 发射度 — 旧口径偏差 +5/+10.5/+25%);
横向 J 为归一化振幅 (× βγ = p_ref/mc, 与 Cemit 文件列单位一致)。
修正后与 golden 对照: z=1.5 位置三平面三分数 **max|dev| = 0.024%**;
z 平面全线 499 位置 max 0.001%; 横向线圈区 (z∈[1.04,1.36] m)
max 0.56% (ASCII 相位 dump 5 位有效数字舍入 + 正则动量重建噪声,
测试按 1% 容差留 2× 余量)。k = round(f·N) 在 N=500 上与
floor/ceil 不可区分 (f·N 恒整数), 保留 round 并注释。
实现: analysis/core_emit.py; 测试: test/test_cemit_cv.py (17 项)。

## 3. slice 发射度与失配参数缺乏 ASTRA 级对照 ★中 (维持)

**现象**: slice 发射度只有解析验证 (高斯束团 slice εn = σxσpx/mc²)
与 ζ≥1 下界; ASTRA 自身不输出 slice 级发射度文件。

**2026-08-19 更新**: 计划的一次性第三方对照 (pmd-beamphysics) 因
.venv 无此包而延后 (R6); 其余链 (统计/slice/核心分数/正电子)
均已有真实运行或解析级验证。若需补: 一次性 pip install +
/tmp 脚本, 不进仓库、不成为运行时依赖。见 codex_tasks.md M7。
**风险**: 低-中 — slice 分析用于展示与失配诊断。

## 4. Error 文件 (ErrorS) 无真实 golden — **已关闭** (2026-08-19, R5)

**结论**: &ERROR (ErrorS=T, Err_MaxB(1)=0.02) 真跑生成 golden
Error.001; read_error 列语义与手册 Table 4 (run#, z, FOM(1..10))
逐列验证; 与名义值运行 (Xemit) 一致性抽查通过。
重要发现: 手册 6.5 的误差按高斯抽样且**无种子参数**, FOM(2) 列
跨 run 不可复现 (三次实测不同), 其余 11 列跨 run 稳定。测试设计:
冻结数值断言仅针对归档快照, 语义断言 (均值≈名义值、散布含名义值、
能量/束长不受螺线管误差影响、FOM(4..10)=0) 对任意新 run 成立。
测试: test/test_error_golden.py (5 项)。

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

## 8. Plasma_Example_2 激光 3D 图头在不同 ASTRA 构建间不可移植 — **已关闭** (2026-08 绝对物理验证通过) ★中

**现象** (2026-08 实测, macOS Apple Silicon 构建): DESY 原版
laser.dat 头三行是 (n, min, spacing) 且计数写成浮点 (8.1e+01):
* 原版浮点计数   -> ASTRA "Error while reading file: laser.dat"
* 逐值网格行 (n 后跟 n 个值) -> 读取 3D 图时 SIGSEGV
* 整数计数 + (n,min,spacing) -> 正常 (fix_laser_map_header 现行
  输出, 已改为此形式)

**归档与运行策略** (R2-3-2): 仓库内 `examples/Plasma_Example_2/laser.dat`
与 DESY 原版逐字节一致 (MD5 68d016175859b20c1e2ccee5057c2d46),
浮点计数头是 DESY 原版特征, 仓库按原版归档; 运行前由
`examples/_examples_spec.py` 的 `laser_fix=True` 做 stage-time 头修复
(`astra_tools.io.field_map.fix_laser_map_header`, 仅头 3 行计数取整,
数据体不动); 整数计数头重跑可字节复现 golden。

**处置**: Plasma_Example_2 的 golden (Xemit/Yemit/Zemit/Log/0011.001)
与 golden_expected.json 已用整数计数头的真实运行重新生成
(末态 E_kin 100 -> 421.8 MeV, eps_nx = 0.156 μm); 旧 golden
(239.6 MeV, eps_nx = 108 μm) 来自旧头形式的运行, 无法在当前构建
复现, 已作废。

**绝对物理验证 (2026-08, 关闭依据)**: 对照 Astra-Manual V3.2 §6.9
线性激光尾场公式独立积分 (n=1e17 cm^-3, λ=800 nm, a0=0.8, σz=6.4 μm,
w0=50 μm):
* 末态动能: 公式积分 422.7 MeV vs golden 421.8 MeV, 吻合 0.2%;
* E_peak = 4.32 GV/m = 0.142·E_WB (线性区, E_WB=30.4 GV/m);
* 失相长度 L_d = λp³/λ² = 1.84 m, 泵浦耗尽 ~1.7 m, 均 >> 0.1 m 靶长;
* 激光功率 ~54 TW << 自聚焦临界功率 ~0.3 PW (P_cr = 17.4·(λp/λ)² GW);
* 束加载效应可忽略 (a0<1 线性区, 尾场对束电荷不敏感);
* √(2π) 系数变体 (末态 743.6 MeV) 与 golden 不符, 被排除。
* 限制: DESY 官方算例说明 PDF 未给出期望末态能量数字, 故以手册
  公式为准 (0.2% 吻合)。详见 docs/dev_manual/codex_tasks.md M6。

## 9. 二进制分布文件无真实 golden — **已关闭** (2026-08-19, R1)

**结论**: Manual_Example 设 BINARY=T 真跑, 归档首个真实二进制分布
文件 golden (examples/Manual_Example/golden/Example_binary.001)。
**格式实测**: 真实 ASTRA 二进制 = Fortran sequential unformatted
记录流 (每条粒子一个 [i32 记录长=72][8×f64 (x,y,z,px,py,pz,clock,
charge) + 2×i32 (species,status)][i32 记录长] 记录; 首条=参考粒子
绝对坐标, 其余 z/pz/clock 相对), **不是**此前假设的 5 值头明文流。
读取器现支持双布局 (真实记录流 + 遗留流), 字节序/种类语义/头 Q=|Q|
修复在新路径上全部生效; 与 ASCII dump 逐粒子对照 + Xemit/Zemit
末行交叉验证通过。测试: test/test_binary_golden.py (5 项)。
注记: 记录流参数 (80/72 字节、i32 记录标记) 已验证于 gfortran
macOS 构建; 其他编译器 (如 ifort 8 字节标记) 回落到遗留路径并
显式报错, 不静默误读 (astra_dist.py 模块 docstring)。

## 建议优先级

剩余开放项: 问题 3 的第三方对照 (R6, 延后, 见 M7)。
其余问题 1/2/4/8/9 已全部关闭; 问题 5 已解释 (非缺陷); 问题 6
(BFF) 无 ASTRA 对照属 ASTRA 不输出 BFF 的固有限制。
详见 docs/dev_manual/codex_tasks.md 第三阶段路线。

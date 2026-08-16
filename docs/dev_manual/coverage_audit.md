# 后处理程序覆盖审计 (手册 5.5-5.7 vs astra-notebook)

> 2026-08 审计。对照 Astra-Manual_V3.2 第 5.5 (lineplot)、5.6 (postpro)、
> 5.7 (fieldplot) 的菜单逐项核对。✅ 已实现; ⚠️ 部分/近似;
> ❌ 未实现 (含理由)。

## lineplot (5.5)

### 菜单 1 (随 z 演化)

| # | 手册项 | 实现 | 状态 |
|---|--------|------|------|
| 1 | 横向发射度 | plot_emittance_evolution | ✅ |
| 2 | rms 束斑 | plot_envelope_evolution | ✅ |
| 3 | rms 发散角 | plot_divergence_evolution | ✅ |
| 4 | 纵向发射度 | plot_lineplot_overview 面板 | ✅ |
| 5 | rms 束长 | plot_bunch_length_evolution | ✅ |
| 6 | rms 能散 | plot_energy_spread_evolution | ✅ |
| 7 | 关联能散 | plot_correlated_energy_spread (新增) | ✅ |
| 8 | 平均能量 | plot_energy_evolution | ✅ |
| 9 | 粒子速度 | plot_velocity_evolution | ✅ |
| 10 | 参考粒子动量 | plot_ref_momentum (新增) | ✅ |
| 11 | dp/dz | plot_ref_trajectory 第 3 面板 | ✅ |
| 12 | 参考粒子轨迹 | plot_ref_trajectory | ✅ |
| 13 | 拉莫尔角 | plot_larmor (avr+rms, 比原程序更全) | ✅ |
| 14/15 | 探针轨迹 (笛卡尔/柱坐标) | plot_probe_trajectories (仅笛卡尔) | ⚠️ |
| 16 | 空间电荷场 | plot_space_charge_fields | ✅ |
| 18 | 含孔径几何 | plot_envelope_with_aperture | ✅ |
| 19 | ZOOM | matplotlib 原生交互缩放 | ✅ |
| 20 | fit/save/read | 导出 CSV/npz 替代; 解析拟合未做 (价值低) | ⚠️ |
| 21/22 | to file / overview | savefig / plot_lineplot_overview | ✅ |
| 23/24 | next/prev run | run_selector 手动选择 (无按钮步进) | ⚠️ |

### 菜单 2 (相位扫描/损失/光学)

| # | 手册项 | 实现 | 状态 |
|---|--------|------|------|
| 1 | Energy vs Phase | plot_phase_scan | ✅ |
| 2 | dE/dz vs Phase | plot_pscan_dedz | ✅ |
| 3 | 压缩因子 (z) | plot_pscan_compression | ✅ |
| 4 | 压缩因子 (时间) | plot_pscan_compression_time (新增) | ✅ |
| 5 | 粒子损失 | plot_losses | ✅ |
| 6 | 能量沉积 | plot_losses (右轴) | ✅ |
| 7 | 束载 | plot_beam_loading | ✅ |
| 9 | 空间电荷缩放因子 | plot_tcheck_scaling | ✅ |
| 10 | 缩放计数器 | plot_tcheck_counter (新增) | ✅ |
| 11 | 平均步长 | plot_step_size_evolution (ref 行距近似) | ⚠️ |
| 12/13 | beta/alpha | plot_beta_alpha | ✅ |
| 14 | 相位推进 | plot_phase_advance | ✅ |
| 15 | 相干长度 | plot_coherence_length | ✅ |

### 菜单 3 (扫描)

| # | 手册项 | 实现 | 状态 |
|---|--------|------|------|
| 1-10 | FOM(1..10) | plot_scan_fom | ✅ |
| 1-10 | Err. FOM 直方图 | plot_error_hist | ✅ |
| 11 | FOM 保存位置 | plot_scan_position (新增) | ✅ |
| 15/16 | bins/title | 函数参数 | ✅ |

### 菜单 4 (缩减/核心/阴极)

| # | 手册项 | 实现 | 状态 |
|---|--------|------|------|
| 1/2 | 缩减发射度 z / z&E | plot_reduced_emittance | ✅ |
| 3 | 发射度差 (标准-缩减) | 未实现 (需同时读 Xemit+Xemit2 对照) | ❌ |
| 4/5 | 相关发射度贡献 (Eq 4.4) | 未实现 (公式拆分未做) | ❌ |
| 6 | 缩减纵向发射度 | 数据列已读 (K2z/K3z), 未画 | ❌ |
| 7/8 | trace-space 发射度 | plot_trace_emittance | ✅ |
| 9-11 | 核心发射度 x/y/z | plot_core_emittance (plane 参数, 新增 y/z) | ✅ |
| 12-14 | 交叉粒子/亚束团 | plot_cr_emit (新增); Sub_emit 读者缺失 | ⚠️ |
| 15/16 | 阴极 Ez / 电荷 | plot_cathode_emission | ✅ |
| 17 | 位置 vs 时间 (头/尾) | 未实现 (Cathode 文件无此列) | ❌ |

## postpro (5.6)

| # | 手册项 | 实现 | 状态 |
|---|--------|------|------|
| 1 | 横向相空间 + 投影 | plot_transverse_phase_space + plot_distributions | ✅ |
| 2 | 纵向相空间 (z) | plot_phase_space(z) | ✅ |
| 3 | 纵向相空间 (时间) | 未实现 (t 坐标版本) | ❌ |
| 4/5 | 三视图 / vs 时间 | plot_overview (x-y/z-x/z-y) | ✅ |
| 6 | 核心发射度 (按粒子数) | plot_core_fraction_curves | ✅ |
| 7 | slice 子菜单 | 见下 | ✅ |
| 8 | 横向核心亮度 | plot_core_brightness | ✅ |
| 9 | 核心束长 | plot_core_fraction_curves (sigma_z 面板) | ✅ |
| 10 | z-plot | plot_z_plot | ✅ |
| 11 | 相空间操作 | cut_distribution (窗口切割); 旋转/鼠标箭头未做 | ⚠️ |
| 14/15 | forward/backward | PhaseStepper 步进器 (03_postpro) | ✅ |
| 17/18 | next/prev run | run_selector (手动) | ⚠️ |

### slice 子菜单 (5.6.3)

| # | 手册项 | 实现 | 状态 |
|---|--------|------|------|
| 1 | slice 发射度 | plot_slice_emittance | ✅ |
| 2 | 失配参数 | plot_slice_mismatch | ✅ |
| 3 | slice 能散 | plot_energy_chirp (平均能量) 近似 | ⚠️ |
| 4 | 投影椭圆 | plot_slice_ellipses_3d (3D) | ⚠️ |
| 5 | z-投影 (px/pz rms vs z) | plot_slice_sizes | ✅ |
| 6 | 3D 椭圆 | plot_slice_ellipses_3d | ✅ |
| 7-9 | 投影切换/按能量切/相关 | 未实现 | ❌ |
| 10 | slice 数 | n_slices 参数 | ✅ |
| 11/12 | 减线性相关/改关联能散 | 未实现 | ❌ |

## fieldplot (5.7)

| # | 手册项 | 实现 | 状态 |
|---|--------|------|------|
| 1 | 腔场 (TM) | plot_cavity_field (Ez/Er/Bphi/半径) | ✅ |
| 2 | TE 模场 | plot_field_profile 通用 2 列 | ⚠️ |
| 3 | 二极模场 | plot_3d_map_slices (TDS 3D 图) | ⚠️ |
| 4 | 螺线管场 | plot_solenoid_field | ✅ |
| 5 | 四极场 | plot_field_profile 通用 | ⚠️ |
| 6/7 | 二极场 (水平/垂直) | plot_field_profile / plot_3d_map_slices | ⚠️ |
| 8 | 弯曲阴极 | plot_curved_cathode_contour | ✅ |
| 9 | 阴极表面场 | plot_cathode_emission (E_spch, 时间曲线) | ⚠️ |
| 10 | 空间电荷场子菜单 | plot_space_charge_fields (简化) | ⚠️ |
| 11 | 含几何 | plot_envelope_with_aperture | ✅ |
| 14 | next page | 多面板静态图 | ✅ |
| (5.7.2) | 3D 场图 | plot_3d_map_slices + plot_laser_on_axis | ✅ |
| (5.7.3) | 激光/等离子体 | plot_laser_on_axis + plot_plasma_profile | ✅ |

## 交互逻辑检查 (2026-08)

* 渲染方式: 2D 相空间为普通散点 (s=8, 确定性子采样 max_points,
  0.5-99.5 百分位裁剪); RMS 椭圆叠加已按用户决定移除。
* 步进: 03_postpro 内置 PhaseStepper; examples/postpro_step_demo.ipynb
  是独立演示 (8 个 z monitor 自动刷新)。
* 功能演示 demo (2026-08 新增, examples/ 下 4 本): generator_demo
  (INPUT 卡改参重跑 + 发射度回读)、bff_demo (直接法/FFT/解析式
  三方对照 + CSR 特征点)、stats_validation_demo (统计量逐列对照
  ASTRA Xemit/Zemit + 螺线管正则动量)、lineplot_demo (lineplot 全
  菜单 + 稀有文件类型合成教学数据)。

* 03_postpro: PhaseStepper 步进 (滑块+◀◀ ◀ ▶ ▶▶) 自动刷新统计与相空间;
  其余单元步进后重跑, 与原 postpro 逻辑一致 ✅
* 04_lineplot: 改为 discover_sim_runs 自动发现 stem (原硬编码 astra,
  换了 deck 名就会 FileNotFoundError — 已修) ✅
* 05_fieldplot: 场文件固定指向 examples/ (教学定位, 注释说明) ✅
* 02_astra 表单 -> 写入: 数组参数数值化、TRUE/FALSE 语义已修 (前轮) ✅
* 空/单粒子/缺文件: 各图有守卫与中文提示 ✅

## 结论

主体功能全覆盖; 剩余 ❌ 项集中在: 解析拟合 (fit/save&read)、
Eq.4.4 相关贡献拆分、纵向时间相空间、slice 减相关/改关联能散、
鼠标交互切割 — 均为低频高级操作, 记录在案, 按需补充。

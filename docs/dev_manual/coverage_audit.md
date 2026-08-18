# 后处理程序覆盖审计 (手册 5.5-5.7 vs astra-notebook)

> 2026-08 审计。对照 Astra-Manual_V3.2 第 5.5 (lineplot)、5.6 (postpro)、
> 5.7 (fieldplot) 的菜单逐项核对。✅ 已实现; ⚠️ 部分/近似;
> ❌ 未实现 (含理由)。

## lineplot (5.5)

> 2026-08-18 复核 + 实施补全: lineplot.ipynb 已接入全部已有函数,
> 后端补纵向 Trace 与 cross-over 束斑面板 (✅=已展示; ❌=未实现)。

### 菜单 1 (随 z 演化)

| # | 手册项 | 实现 | 状态 |
|---|--------|------|------|
| 1 | 横向发射度 | plot_emittance_evolution | ✅ |
| 2 | rms 束斑 | plot_envelope_evolution | ✅ |
| 3 | rms 发散角 | plot_divergence_evolution (overview 面板) | ✅ |
| 4 | 纵向发射度 | plot_lineplot_overview 面板 | ✅ |
| 5 | rms 束长 | plot_bunch_length_evolution (t 版已调, z 版在 overview) | ✅ |
| 6 | rms 能散 | plot_energy_spread_evolution (overview 面板) | ✅ |
| 7 | 关联能散 | plot_correlated_energy_spread | ✅ |
| 8 | 平均能量 | plot_energy_evolution | ✅ |
| 9 | 粒子速度 | plot_velocity_evolution | ✅ |
| 10 | 参考粒子动量 | plot_ref_momentum | ✅ |
| 11 | dp/dz | plot_ref_trajectory 第 3 面板 | ✅ |
| 12 | 参考粒子轨迹 | plot_ref_trajectory | ✅ |
| 13 | 拉莫尔角 | plot_larmor (avr+rms) | ✅ |
| 14 | 探针轨迹 (笛卡尔) | plot_probe_trajectories | ✅ |
| 15 | 探针轨迹 (柱坐标 r/z, x/y) | plot_probe_trajectories(mode='cylindrical') (r(z) + x/y 投影) | ✅ |
| 16 | 空间电荷场 | plot_space_charge_fields | ✅ |
| 18 | 含孔径几何 | plot_envelope_with_aperture | ✅ |
| 19 | ZOOM | matplotlib 原生交互缩放 | ✅ |
| 20 | fit, save & read | 导出 CSV/npz 替代; 解析拟合未做 (价值低) | ⚠️ |
| 21/22 | to file / overview | savefig / plot_lineplot_overview | ✅ |
| 23/24 | next/prev run | run_selector 手动选择 (无按钮步进) | ⚠️ |

### 菜单 2 (腔相位/损失/空间电荷/光学函数)

| # | 手册项 | 实现 | 状态 |
|---|--------|------|------|
| 1 | Energy vs. Phase | plot_phase_scan | ✅ |
| 2 | dE/dz vs. Phase | plot_pscan_dedz | ✅ |
| 3 | Compression factor (z) | plot_pscan_compression | ✅ |
| 4 | Compression factor (time) | plot_pscan_compression_time | ✅ |
| 5 | Particle loss | plot_losses (lost/m) | ✅ |
| 6 | Energy deposition | plot_losses 双轴 (J/m) | ✅ |
| 7 | Beam loading | plot_beam_loading | ✅ |
| 9 | Sp. ch. scaling factors | plot_tcheck_scaling | ✅ |
| 10 | Sp. ch. scaling counter | plot_tcheck_counter | ✅ |
| 11 | Average time step | plot_step_size_evolution | ✅ |
| 12/13 | β/α 函数 | plot_beta_alpha | ✅ |
| 14 | phase advance | plot_phase_advance | ✅ |
| 15 | coherence length | plot_coherence_length | ✅ |

### 菜单 3 (扫描 FOM / 误差直方图)

| # | 手册项 | 实现 | 状态 |
|---|--------|------|------|
| 1-10 左 | FOM(1-10) | plot_scan_fom (i 参数逐个) | ✅ |
| 1-10 右 | Err. FOM 直方图 | plot_error_hist | ✅ |
| 11 左 | FOM 保存位置 | plot_scan_position | ✅ |
| 11 右 | Err. 保存位置 | 未实现 | ❌ |
| 15/16 | bins / title | GUI 交互 (n_bins 参数) | ⚠️ |

### 菜单 4 (缩减/核心/Trace/cross-over 发射度)

| # | 手册项 | 实现 | 状态 |
|---|--------|------|------|
| 1 | reduced trans. emittance z | plot_reduced_emittance (eps_red_z) | ✅ |
| 2 | reduced trans. emittance z&E | plot_reduced_emittance (eps_red_zE) | ✅ |
| 3 | trans. emittance difference | plot_emittance_difference (标准-缩减, z 插值对齐; 无 Yemit2 时 y 用 x 近似标注 est.) | ✅ |
| 4/5 | cor. emittance contributions x/y | plot_correlated_emittance_contributions (Eq. 4.4 K2z/K3z/K2E/K3E) | ✅ |
| 6 | reduced long. emittance | plot_reduced_longitudinal_emittance (当前分布, 减去 2nd/3rd 阶 z-pz 相关, 手册 4.13.6) | ✅ |
| 7 | trans. Trace Space emittance | plot_trace_emittance (eps_tr_x/y) | ✅ |
| 8 | long. Trace Space emittance | plot_trace_emittance 纵向 eps_tr_z (2026-08 新增) | ✅ |
| 9-11 | core emittance x/y/z | plot_core_emittance (C95/C90/C80) | ✅ |
| 12 | emittance w.o. cross over | plot_cr_emit (eps_x/y + 电荷) | ✅ |
| 13 | beam size w.o. cross over | plot_cr_emit 束斑面板 x_rms/y_rms (2026-08 新增) | ✅ |

## postpro (5.6)

> 2026-08-18 全覆盖实施完成 (见 postpro_full_coverage_plan.md)。notebook
> postpro 已接入全部新增功能。

| # | 手册项 | 实现 | 状态 |
|---|--------|------|------|
| 1 | 横向相空间 + 投影 | plot_transverse_phase_space (x-x'/y-y') + plot_distributions (x/y/z 投影) | ✅ |
| 2 | 纵向相空间 (z) | plot_phase_space(z) | ✅ |
| 3 | 纵向相空间 (时间) | plot_phase_space(plane="t") (bunch_time, 手册 5.6 时间坐标规则) | ✅ |
| 4 | Front/Top/Side 三视图 | plot_overview (3x2 超集) | ✅ |
| 5 | 三视图 vs 时间 | plot_overview(time=True) | ✅ |
| 6 | 核心发射度 (按粒子数) | compute_core_emittance_by_fraction + plot_core_emittance_curve (4.13.5 单粒子振幅排序; f=1 与 ASTRA 精确一致, f<1 趋势一致) | ✅ |
| 7 | slice 子菜单 | 见下 | ✅ |
| 8 | 横向核心亮度 | plot_core_brightness (postpro 已接入, 需 Xemit/Yemit/LandF) | ✅ |
| 9 | 核心束长 | plot_central_charge_fraction_curves (sigma_z 面板) | ✅ |
| 10 | z-plot (含丢失粒子) | plot_z_plot (active/lost/passive 分色) | ✅ |
| 11 | 相空间操作 | optimized_cut + CutControls (滑块+应用/撤销) + rotate_phase_space + cut_distribution | ✅ |
| 13 | fit, save & read | OverlayManager 内存叠加 (无参数拟合) | ⚠️ |
| 14/15 | forward/backward | PhaseStepper 步进器 (postpro) | ✅ |
| 16 | to file | export_distribution / export_statistics | ⚠️ |
| 17/18 | next/prev run | 手动选择 | ⚠️ |

### 菜单 2 (5.6.2 任意参数对)

| # | 手册项 | 实现 | 状态 |
|---|--------|------|------|
| - | 任意 abscissa/ordinate 参数选择 | param_columns (12 参数) + plot_arbitrary + PhaseSpaceParamSelector | ✅ |
| - | add projections | plot_arbitrary(add_proj=True) (边缘直方图) | ✅ |
| - | subtract linear correlation | plot_arbitrary(subtract_corr=True) | ✅ |
| - | 叠加图 (save & read overlay) | OverlayManager (plot_arbitrary 模块) | ✅ |

### slice 子菜单 (5.6.3)

| # | 手册项 | 实现 | 状态 |
|---|--------|------|------|
| 1 | slice 发射度 | plot_slice_dashboard / plot_slice_emittance | ✅ |
| 2 | 失配参数 | plot_slice_mismatch | ✅ |
| 3 | slice 能散 | plot_energy_chirp (平均能量 + sigma_E/E) | ✅ |
| 4 | 投影椭圆 | plot_slice_ellipses_2d (2D 投影 rms 椭圆; 项 7 平面切换 xxp/yyp/xyp/yxp, 项 11 减线性相关) | ✅ |
| 5 | z-投影 (x_rms + 发散角 rms/avr vs z) | plot_slice_sizes(divergences=True) (含 σx'/σy' 双轴) | ✅ |
| 6 | 3D 椭圆 | plot_slice_ellipses_3d | ✅ |
| 7 | 投影切换 (x-px/y-py/x-py/y-px) | plot_slice_ellipses_3d(plane="xxp"/"yyp"/"xyp"/"yxp") | ✅ |
| 8/9 | w.r.t. position / w.r.t. energy 切片 | compute_slice_analysis(binning="equi_energy") | ✅ |
| 10 | slice 数 | n_slices 参数 (equi_spaced / equi_charge / equi_energy) | ✅ |
| 11 | 减线性相关 (slice 椭圆) | plot_slice_ellipses_3d(subtract_corr=True) | ✅ |
| 12 | 改关联能散 | modify_correlated_energy_spread(dist, factor) | ✅ |
| 15/16 | 数据/图 to file | export_* 近似 | ⚠️ |

### 相空间操作 (5.6.4) 与通用规则

| 手册项 | 实现 | 状态 |
|--------|------|------|
| 旋转到旋转坐标系 | rotate_phase_space (x-y 平面旋转) | ⚠️ |
| 优化切割 (interval 参数) | optimized_cut / optimized_cut_center (存活粒子最大化) | ✅ |
| 鼠标箭头切割 | CutControls (滑块替代鼠标, 应用/撤销) | ✅ |
| accept/reject changes | CutControls 撤销 (保留原分布副本) | ✅ |
| 保存新分布供继续追踪 | write_distribution (二进制/ASCII, 与 reader 互逆) | ✅ |
| 状态过滤规则 (≥-6 绘图, z-plot 全粒子, slice/core 用 status>1) | color_by_status 用 ≥-6; active=status>1; z-plot 含 lost/passive | ✅ |
| Table 6 状态颜色/符号编码 | scatter2d(color_by_status=True) (8 类颜色/符号, 图例) | ✅ |
| Plot_steering.par (重定向统计/显示条件) | read_plot_steering (Stat/ion_mass/CP_ind 解析) | ✅ |
| 混合分布 CP_ind 粒子索引着色 | cp_index_colors (颜色数组, 供绘图) | ✅ |
| 显示 generator 的 .ini 分布 | read_distribution 支持 | ✅ |

## fieldplot (5.7)

> 2026-08-18 复核 + 实施补全: fieldplot 已接入全部已有函数 (阴极/激光
> 激活, 弯曲阴极/场剖面/空间电荷新增), 均带文件存在守卫。✅=已展示;
> ⚠️=简化/部分; ❌=未实现。2026-08 (遗留项收尾): 新增 TE 模场 (Bz/Br/Eφ),
> 场展开半径 R3rd (TM/TE/螺线管, 手册 8 章), 螺线管/四极 next page,
> 激光 rms 包络+焦点, 弯曲阴极电荷环, 等离子体场 vs z/zeta。

| # | 手册项 | 实现 | 状态 |
|---|--------|------|------|
| 1 | 腔场 (TM): Ez/Er/Bφ 振幅 | plot_cavity_field (2D 图 + 轴上) | ✅ |
| 1' | 腔场 next page: 场展开半径 R3rd | plot_field_expansion_radius (TM: R3rd_Er/Bφ, 手册 8 章; 数值噪声处 NaN 跳过) | ✅ |
| 2 | TE 模场: Bz/Br/Eφ 振幅 | TEField + read_te_field + plot_te_field (Bz/Br/Eφ 2D + 轴上 Bz; 手册 6.9 TE_ 前缀 + 8 章展开) | ✅ |
| 3 | 二极模场 (TDS) | plot_3d_field_map (3D 图) | ⚠️ |
| 4 | 螺线管场: Bz + 径向梯度 | plot_solenoid_field (Bz 2D + Br 箭头) | ✅ |
| 4' | 螺线管 next page: 单独 Br / R3rd | plot_solenoid_components (单独 Bz/Br + R3rd) + plot_field_expansion_radius (静磁) | ✅ |
| 5 | 四极场 (水平/垂直梯度) | plot_field_profile 通用 1D | ✅ |
| 5' | 四极 next page: Bz + 综合图 | plot_quadrupole_field (Gx/Gy 主图 + Bz/综合 next page; 理想四极 Bz=0 标注) | ✅ |
| 6/7 | 二极场 (水平/垂直) | plot_field_profile / plot_3d_field_map | ✅ |
| 8 | 弯曲阴极轮廓 | plot_curved_cathode_contour (列解读已修复, 手册 4.4.5) | ✅ |
| 8' | 弯曲阴极: 电荷环位置 | plot_curved_cathode_contour(show_rings=True) (环在表面背面, 手册 4.4.5) | ✅ |
| 9 | 阴极表面场 | plot_cathode_emission (E_acc/E_spch/电荷 vs t; 无 Cathode 文件走合成) | ✅ |
| 10 | 空间电荷场子菜单 (Er/Ez/Bfi/Er eff/电荷密度/线电荷密度/速度剖面) | plot_space_charge_fields (仅探针 Er/Ez 简化) | ⚠️ |
| 11 | 含几何 | plot_envelope_with_aperture | ✅ |
| 14 | next page | 多面板静态图 | ✅ |
| (5.7.2) | 3D 场图 (1-6) | plot_3d_field_map (矢量剖面/2D 等值线/3D 等值线/标量剖面) | ✅ |
| (5.7.3) | 激光 3D 图截面 (1-6) | plot_laser_on_axis (轴上剖面 vs z/t; 无真实文件走合成) | ✅ |
| (5.7.3)' | 激光: rms 束包络 + 焦点位置 | plot_laser_envelope (f² 加权 rms 包络 σx/σy + 焦点标注, 5.7.3) | ✅ |
| (5.7.3) | 等离子体 fields vs z (项 8) | plot_plasma_fields(vs='z') (线性尾场解析模型 + 密度, 手册附录) | ✅ |
| (5.7.3) | 等离子体 fields vs zeta (项 9) | plot_plasma_fields(vs='zeta') (共动参数 ζ = z - c t) | ✅ |

## 交互逻辑检查 (2026-08)

* 渲染方式: 2D 相空间为普通散点 (s=8, 确定性子采样 max_points,
  0.5-99.5 百分位裁剪); RMS 椭圆叠加已按用户决定移除。
* 步进: postpro 内置 PhaseStepper; examples/postpro_demo.ipynb
  是独立演示 (8 个 z monitor 自动刷新)。
* 功能演示 demo (2026-08 新增, examples/ 下 6 本): generator_demo
  (INPUT 卡改参重跑 + 发射度回读)、bff_demo (直接法/FFT/解析式
  三方对照 + CSR 特征点)、stats_validation_demo (统计量逐列对照
  ASTRA Xemit/Zemit + 螺线管正则动量)、lineplot_demo (lineplot 全
  菜单 + 稀有文件类型合成示例数据)、fieldplot_demo (腔场二维图 +
  3D 场图四种视图 + 螺线管 + 阴极/激光/等离子体场剖面)、
  postpro_demo (postpro 全菜单 + 步进演示)。

* postpro: PhaseStepper 步进 (滑块+◀◀ ◀ ▶ ▶▶) 自动刷新统计与相空间;
  其余单元步进后重跑, 与原 postpro 逻辑一致 ✅
* lineplot: 改为 discover_sim_runs 自动发现 stem (原硬编码 astra,
  换了 deck 名就会 FileNotFoundError — 已修) ✅
* fieldplot: 场文件固定指向 examples/ (示例定位, 注释说明) ✅
* astra 表单 -> 写入: 数组参数数值化、TRUE/FALSE 语义已修 (前轮) ✅
* 空/单粒子/缺文件: 各图有守卫与中文提示 ✅

## 结论

主体功能全覆盖; 剩余 ❌ 项集中在: 解析拟合 (fit/save&read)、
Eq.4.4 相关贡献拆分、纵向时间相空间、slice 减相关/改关联能散、
鼠标交互切割 — 均为低频高级操作, 记录在案, 按需补充。

# postpro 全覆盖方案 (postpro vs 手册 5.6)

> 2026-08-18 规划。目标: 让 postpro 覆盖原程序 postpro (手册 5.6) 的
> 全部功能。后端逻辑在 `astra_tools/`, notebook 只编排; 遵循 AGENTS.md
> 物理规则与五层测试纪律。交互式功能 (鼠标切割/叠加) 用 ipywidgets
> 简化, 保持可用即可。

## 0. 现状摘要 (已覆盖)

横向/纵向相空间(z)、三视图(超集)、z-plot(含丢失)、slice 发射度/失配/
能散/3D 椭圆/切片数、PhaseStepper 步进、导出 —— 均 ✅。
本方案只补 ❌/⚠️ 项。

## 1. 时间坐标视图 (菜单1 项3 纵向相空间·时间, 项5 三视图 vs 时间)

- 新增 `astra_tools/analysis/time.py`:
  - `bunch_time(dist)` -> np.ndarray [s]: 手册 5.6 开头规则 ——
    not_started (status -1..-6) 用发射时间 `clock`; 其余用
    t = (z - <z>) / (β̄·c), β̄ 由平均纵向动量算; 混合分布发 warning。
  - `bunch_time_ps(dist)` 便捷 ns/ps 版。
- `plot_phase_space(dist, plane="t")`: 纵向相空间时间版 (t vs dp/p)。
- `plot_overview(dist, time=True)`: 三视图 vs 时间 (t-x / t-y / t-x')。
- postpro 相空间 cell 加 "z / t" 切换。

## 2. 菜单2: 任意参数对相空间 (5.6.2)

- 新增 `astra_tools/plot/arbitrary_phase_space.py`:
  - `param_columns(dist)` -> {名称: (数组[SI], 单位串)}: x, y, z, px, py,
    pz, clock, x', y', dp/p, E_kin, t。
  - `plot_arbitrary(dist, xp, yp, subtract_corr=False, add_proj=False,
                    color_by_status=False)`: 任意两参数散点;
    subtract_corr = 拟合 v = a + b·u 后画残差 (叠加线性相关去掉)。
- `astra_tools/widgets/selectors.py` 新增 `PhaseSpaceParamSelector`:
  两个参数下拉 + "加投影" "减线性相关" "状态着色" 复选框 + 叠加。
- 叠加 (overlay): `PlotOverlayManager` 内存保存当前散点数据,
  叠加时按 Plot_mode 颜色反转或黑色。

## 3. slice 增强 (5.6.3)

- 项 8/9: `compute_slice_analysis(binning="equi_energy")` —— 按能量
  等宽切片 (复用 equi_charge 的边沿修复逻辑)。
- 项 5: `compute_slice_analysis` 补充 slice 发散角统计
  (x'_rms, y'_rms, x'_avr, y'_avr, 正则动量/p_ref); `plot_slice_sizes`
  增加 x'/y' 曲线。
- 项 7: `plot_slice_ellipses_3d(plane=...)` 支持 (x,x')/(y,y')/(x,y')/
  (y,x') 四种投影切换。
- 项 11: `plot_slice_ellipses_3d(subtract_corr=True)` 椭圆减线性相关。
- 项 12: `modify_correlated_energy_spread(dist, factor)` 缩放 pz 的
  关联 (随 z) 部分, 重算 slice 观察变化。

## 4. 相空间操作 (5.6.4)

- `astra_tools/analysis/cuts.py`:
  - `optimized_cut(dist, param, target)`: 优化切割 —— 给定 interval 参数
    (如 z 束长/x 尺寸), 找使存活粒子数最大的对称窗口。
  - 旋转: `rotate_phase_space` 已有 (x-y 平面)。
- `astra_tools/io/astra_dist.py`: 新增 `write_distribution(dist, path,
  format="ascii")` 保存新分布供继续追踪 (表头: N 粒子 + 参考粒子绝对
  坐标, 其余相对坐标, 与 io/astra_dist.py 读取约定互逆)。
- `astra_tools/widgets/panels.py`: `CutControls` (x/y/z/E 滑块 + 应用/
  撤销, 撤销 = 保留原分布副本)。
- postpro 新增 cell。

## 5. 核心发射度 Cemit 口径 (菜单1 项6, 手册 4.13.5)

> 现状: `read_cemit_file` 已有 (ASTRA 输出, 已金样验证); `plot_core_emittance`
> 是 lineplot 用 (eps_n + C95/90/80 vs z); `plot_central_charge_fraction_curves` 是
> 中心电荷分数口径, **非** 4.13.5。

- 新增 `astra_tools/analysis/core_emit.py`:
  - `single_particle_emittance_contributions(dist, plane)`: 各粒子对
    rms 发射度的单粒子贡献 ε_i (4.13.5, 由 Σ 矩阵分解)。
  - `compute_core_emittance_by_fraction(dist, fractions=(0.5..1.0))`:
    贡献升序排序后按百分比截取, 返回各分数下的 eps_n (x/y/z) 累计曲线。
- `plot_core_emittance_curve(dist)`: 发射度 vs 粒子百分比 (P0)。
- 验证: 与 `read_cemit_file` 的 C95/C90/C80 金样对比, 0.5% 内。
- postpro: 新增 cell, 与现有中心电荷分数图并列 (标注口径差异)。

## 6. 核心亮度接入 notebook (菜单1 项8)

- `plot_core_brightness(ce, landf)` 已存在; 检查其输入 (需 emit 对象 +
  LandF), 在 postpro 新增 cell 用当前 dist + 读 LandF 调用。
- 若无 LandF, 用归一化亮度降级 (函数已处理)。

## 7. 状态着色 Table 6 (通用)

- `astra_tools/plot/phase_space.py::scatter2d` 增加
  `color_by_status=False` 参数: 按手册 Table 6 八类分组着色:
  secondary(>5) 绿 / normal(2,3,5) 黑 / marked(4) 红 / passive(0,1)
  绿+ / cathode(-1..-6) 棕 / lost-aperture(-12..-25) 红点 /
  lost(-26..-30) 蓝圈 / lost(≤-30) 蓝星。
  布局限制用颜色 + 图例 (不做 PGPLOT 符号, 用 marker 区分部分)。
- 默认关闭, 在任意参数对面板与相空间图可开启。

## 8. Plot_steering.par + CP_ind (通用, P2)

- 新增 `astra_tools/io/plot_steering.py`:
  - `read_plot_steering(path)`: 解析 namelist `Steering_parameters`
    (Stat(2,-100:100) 重定向、ion_mass、CP_ind_1..15 RGB)。
  - 将 Stat 重定向应用到统计/显示条件 (仅 slice/core/相空间切割)。
- CP_ind 混合分布按粒子索引着色 (配合 Plot_mode=1)。

## 9. 优先级与里程碑

| 里程碑 | 内容 | 工作量 |
|--------|------|--------|
| M1 (P0) | 时间坐标(项3/5) + Cemit 累计曲线(项6) + 状态着色(Table 6) | 中 |
| M2 (P0/P1) | 核心亮度接入(项8) + slice 能量切片(项8/9) + 发散角 z-投影(项5) | 中 |
| M3 (P1) | 任意参数对(5.6.2) + 叠加 + slice 投影切换/减相关/改关联能散(项7/11/12) | 大 |
| M4 (P1) | 相空间操作: 优化切割/保存分布/CutControls (5.6.4) | 大 |
| M5 (P2) | Plot_steering.par + CP_ind | 小 |

## 10. 文件清单

新增:
- `astra_tools/analysis/time.py`
- `astra_tools/analysis/core_emit.py`
- `astra_tools/plot/arbitrary_phase_space.py`
- `astra_tools/io/plot_steering.py`
- `test/test_core_emit.py`, `test/test_time_coords.py`, `test/test_cuts_optimized.py`

修改:
- `astra_tools/analysis/slices.py` (能量切片 + 发散角)
- `astra_tools/plot/slice_plots.py` (项5/7/11)
- `astra_tools/plot/phase_space.py` (t 平面 + 状态着色)
- `astra_tools/plot/overview.py` (time=True)
- `astra_tools/plot/advanced_plots.py` (Cemit 曲线 / 亮度接入辅助)
- `astra_tools/analysis/cuts.py` (optimized_cut)
- `astra_tools/io/astra_dist.py` (write_distribution)
- `astra_tools/widgets/selectors.py` (PhaseSpaceParamSelector / CutControls)
- `astra_tools/export.py` (可选: 保存新分布入口)
- `notebooks/postpro.ipynb` (接入全部新功能)
- `docs/dev_manual/coverage_audit.md` (完成时回填状态)

## 11. 测试计划 (五层)

1. 单测: time 坐标计算 / Cemit 单粒子贡献 / 能量切片边界 / 优化切割
   存活粒子最大化 / write-read 往返。
2. 金样回归: compute_core_emittance_by_fraction vs `read_cemit_file`
   C95/90/80 (0.5%)。
3. 交叉验证: slice 发散角 vs Xemit x'_rms。
4. 绘图断言: 每个新图轴标签/单位 (test_plots.py 扩展)。
5. e2e: nbconvert 执行 postpro (需 ASTRA/Generator 二进制)。

## 12. 风险与取舍

- 鼠标切割用滑块替代 (ipywidgets), 不做真鼠标交互 (notebook 自然)。
- 单粒子发射度贡献分解对混合电荷符号需用 |q| 权重 (AGENTS.md 规则 7)。
- Cemit 口径与本仓库"中心电荷分数"并存, 标注差异, 不删除旧图。
- 状态着色符号受 matplotlib 限制, 以颜色 + 图例为主, 符号为辅。

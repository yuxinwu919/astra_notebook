# 07 能散: σ_E/E 从全动量动能计算

## 公式
每粒子 E_kin = √((px²+py²+pz²) + m²c⁴) − mc² (**全动量**, 2026-08 审计
P2-2: 旧 pz-only 口径对大发散束低估 mean_E 可达 ~60%);
σ_E/E = std(E_kin)/mean(E_kin)  (不是 σ_p/p)。
σ_p/p 单独保留 (两者在低能区相差 β² 因子)。

## 与 ASTRA 口径的确认 (2026-08 第三阶段 R2b, 决定性)
90deg 大发散探针束真跑 (px≈3.5 MeV/c, pz≈1.87 MeV/c, |p|=3.975 MeV/c)
对照 ASTRA Zemit:
* mean_E:  全动量 3.4967 MeV = Zemit 3.4967 (pz-only 1.4256, 差 59%);
* σ_E:     全动量 21.180 keV = Zemit 21.180 (pz-only 27.763, 差 31%);
* eps_z:   全动量 94.872 eV·m = Zemit 0.094872 (pz-only 差 33%);
* 第 7 列协方差 <zE'>_avr = 21.179 keV 亦为全动量值。
结论: **ASTRA Zemit 的每粒子动能采用全动量 |p| 口径**, 与本实现一致
(test/test_90deg_energy_cv.py 固化, 判别力断言防退化)。

## 低能区验证 (2026-08 R2a)
5 MeV (γ=10.78) 与 5.1 keV (γ=1.00998) 双束真跑, 全列对照 < 4e-4;
γ≈1.01 束上 βγ 归一化发射度 eps_nx = 0.136166 vs ASTRA 0.136160 —
"γ 从动量算"口径在低能区成立 (若误把 p 当动能, βγ 差 3.9 倍,
测试判别力断言排除之)。注意: 低能束 px/p0≈0.2%, 不构成 σE 口径的
独立判别 (该判别完全依赖 90deg 探针束)。
已知细微差异: lowg 束 sig_xp 偏差 4.1e-4 (容差 5e-4, 裕度 1.2×) —
ASTRA 内部发散度定义在低能区与本实现 (px/p_ref) 有 ~4e-4 级差异,
ASTRA 值居中, 与既有 test_cross_validation.py 同源容差一致。

## 数值验证 (1 GeV, 低发散)
E_kin = 999.989145 MeV vs Zemit 999.990000;
σ_E = 1.466645 keV vs Zemit 1.466700。

## 代码位置
astra_tools/constants.py (kinetic_energy_from_momentum_vector);
astra_tools/analysis/statistics.py 及 slices/cuts/core_emit/绘图入口
统一使用全动量口径。

# 07 能散: σ_E/E 从动能计算

## 公式
每粒子 E_kin = √(pz² + m²c⁴) − mc²;
σ_E/E = std(E_kin)/mean(E_kin)  (不是 σ_p/p)。
σ_p/p 单独保留 (两者在低能区相差 β² 因子)。

## 数值验证
E_kin = 999.989145 MeV vs Zemit 999.990000;
σ_E = 1.466645 keV vs Zemit 1.466700。

## 代码位置
astra_tools/analysis/statistics.py (kinetic_energy_from_momentum)

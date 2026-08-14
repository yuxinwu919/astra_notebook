# 04 螺线管中的正则动量发射度

## 公式 (手册 §4.13.1)
p̃x = px + e·Bz·y/2 → [eV/c]: p̃x = px + c·Bz·y/2
p̃y = py − e·Bz·x/2 → [eV/c]: p̃y = py − c·Bz·x/2
Bz = 束团中心处轴上螺线管场; 散角 x' = p̃x / p_ref。

## 数值验证 (Manual_Example, 螺线管 MaxB=0.35 T @ z=1.2 m, 束团 z=1.5 m)
- trace-space ε_geom = 5.32e-3 m·rad (错 7 个量级)
- 正则动量 ε_n = 1.000128 mm·mrad vs Xemit 1.000300 (0.017%)
- σ_x' = 0.743526 vs Xemit 0.743650 mrad

## 代码位置
astra_tools/analysis/emittance.py (canonical_divergence),
statistics.py (bz_on_axis_T 参数), test/test_cross_validation.py

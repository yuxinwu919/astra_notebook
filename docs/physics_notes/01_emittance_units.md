# 01 发射度单位: π mm mrad ≡ mm mrad

## 公式
RMS 相空间椭圆 xᵀΣ⁻¹x = 1 的半轴 a,b 满足 a·b = ε_rms = √det Σ;
椭圆面积 = π·a·b = π·ε_rms。ASTRA 输出的「发射度」是椭圆面积,
故单位写作 π mm mrad, π 被吸收进单位名。

**数值规则**: 文件值 × 1e-6 = ε_rms [m·rad]。绝不额外乘 π。
纵向: π keV mm, 数值 = keV·mm = eV·m (×1)。

## 手册出处
ASTRA Manual V3.2 §4.13、Table 4 (Xemit/Zemit/Cemit 单位栏)。

## 数值验证
- Generator Nemit_x=1 → 实测 σ_px = 681.4 eV/c → ε_n = 1.000009e-6 m·rad = 1.000009 mm·mrad
- Xemit 第 6 列 = 1.0003 (z=1.5 m) = 我的独立计算 ×1e6 (0.017%)
- Nemit_x=2 → 列值 2.0000→1.9993, 线性无 π 因子
- 无螺线管算例: Xemit 列 1.0006 = 独立计算 1.000634 mm·mrad (0.003%)
- lume-astra parser 单位表: 'mm-mrad' 因子 1e-6; pmd-beamphysics: √det(cov)/mc² 单位米

## 代码位置
astra_tools/constants.py (EMIT_M_TO_MM_MRAD=1e6), io/astra_emit.py,
plot/emit_plots.py (轴标签 [π mm mrad]), test/test_cross_validation.py

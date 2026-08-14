# 06 Sigma 矩阵归一化之谜 (待定)

## 已知事实
Sigma 文件 23 列: z[m], E_kin[MeV], 6x6 协方差上三角 21 元素
(手册 Table 4)。坐标系: (x, p̃x, y, p̃y, z, E_kin) (手册 §4.13)。

## 实证
- sig(1,1) = σ_x² [m²] 与 Xemit 精确一致
- sig(5,5) = σ_z² [m²] 与 Zemit 精确一致
- sig(2,2) 与 sig(6,6) 比 σ_x'² / σ_E² 大 ~3.83 倍, 原因不明
  (可能: 未中心化、γ 加权、或坐标定义与手册行文不同)

## 处理
本库的 eigen-emittance 标记为实验性; 验证过的发射度一律用 Xemit。
待数值实验 (对比多个 z 位置、不同束流) 解开因子。

## 代码位置
astra_tools/io/astra_emit.py (read_sigma_file 注释)

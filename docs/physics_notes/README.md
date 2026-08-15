# 物理审查备忘录

每条约定按「公式 → 手册出处 → 数值验证」三段记录。全部验证基于仓库内
真实 ASTRA 输出 (examples/Manual_Example 与 data/emit_unit_check/run1-3),
与 ASTRA 自身的 Xemit/Zemit 输出直接比对, 容差 < 0.5%。

| # | 条目 | 状态 |
|---|------|------|
| 01 | 发射度单位: π mm mrad ≡ mm mrad | 已定案 |
| 02 | active 粒子: status > 1 | 已定案 |
| 03 | 分布文件坐标: z/pz/clock 相对参考粒子 | 已定案 |
| 04 | 螺线管中的正则动量发射度 | 已定案 |
| 05 | Xemit 相关列 = cov(u,u')/σ_u | 已定案 |
| 06 | Sigma 矩阵归一化 (动量/mc, 能量/mc^2; 3.83=1/mc^2) | 已定案 |
| 07 | 能散 σ_E/E 从动能计算 | 已定案 |
| 08 | Zemit 相关列单位 (keV, cov(z,E)/σz) | 已定案 |
| 09 | PScan/Scan 真实数据交叉验证 | 已定案 |
| 10 | 未解决的物理问题报告 (Sigma 3.83 / Cemit / Error / Plasma_2 激光等) | 开放 |
| 11 | IO 格式要点 (3D 场图双头/尾场多块/track 8 列/9-10 列歧义) | 已定案 |

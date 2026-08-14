# 03 分布文件坐标: z/pz/clock 相对参考粒子

## 公式
手册 Table 1: 「Longitudinal particle coordinates, i.e. z, pz and t are
given relative to the reference particle」。ASCII 文件首行 = 参考粒子
(绝对坐标), 其余行 z/pz/clock 为相对值; 二进制文件以 5 值头提供
t_ref/p_ref, 粒子 z/pz/clock 仍为相对值。读取时统一转换为绝对坐标。

## 数值验证
- 不归一化: 499/500 粒子 pz≈±2 keV/c、σ_E=44.7 MeV (荒谬)
- 归一化后: σ_E=1.466645 keV 与 Zemit 1.466700 吻合 0.004%
- σ_z=0.659523 mm 与 Zemit 0.659520 吻合

## 代码位置
astra_tools/io/astra_dist.py (ASCII 首行参考粒子 + 相对→绝对)

# 05 Xemit 第 7 列 = cov(u,u')/σ_u

## 公式
列名 (lume-astra): cov_x__xp/sigma_x, 单位 mrad, 因子 1e-3。
cov(u,u') = 列值 × 1e-3 × σ_u。

## 数值验证 (无螺线管算例)
列值 −3.9889e-5 mrad × σ_x(7.471e-4 m) = −2.980e-11 m·rad
= 独立计算 cov(x,x') ✓ 精确吻合。

## 代码位置
astra_tools/io/astra_emit.py (parse_output_file standardize_labels)

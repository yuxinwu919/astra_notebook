# 02 active 粒子判据: status > 1

## 公式
手册 §4.13: 「All particles with status flag > 1 are taken into account
for the emittance calculation」。

| status | 含义 |
|--------|------|
| > 1 | active (统计/发射度计入; 含轨迹探针 3、标准粒子 5、二次电子) |
| 0, 1 | passive (被追踪但统计排除) |
| -6..-1 | 阴极未发射 |
| < -6 | 丢失 |

旧代码 (beamscope/utils) 用 status == 0 判活跃 — 真实束团 status=5
会被判为全丢失, 是灾难性错误。

## 数值验证
Example.0150.001: 494×status5 + 6×status3; 以 status>1 计算 σ_x =
0.7499624 mm 与 Xemit 0.7499600 mm 吻合 0.0003%。

## 代码位置
astra_tools/distribution.py (active/passive/not_started/lost)

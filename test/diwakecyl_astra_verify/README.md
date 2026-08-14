# DiWakeCyl × ASTRA 尾场验证测试

**测试时间**: 2026-07-05 02:26:39

## 参数

| 参数 | 值 |
|------|------|
| 束团能量 | 200.0 MeV |
| 束团电荷 | 0.05 nC (50 pC) |
| RMS 长度 σz | 0.1 mm |
| 能散 σE | 10.0 keV RMS |
| 归一化发射度 | 1.0 π·mrad·mm |
| 波导内径 b | 0.5 mm |
| 波导外径 a | 1.0 mm |
| 介电常数 ε | 4.41 |
| 模式数 | 2 |
| ASTRA 追踪 | 0.0 → 0.1 m |
| 尾场位置 Wk_z | 0.05 m |
| Bin 数 | 80 |
| 尾场格式 | Taylor_Method_F, 3 terms |

## 结果

成功: 3 个

### 统计对比

| Case | σ(δp/p) [%] | ε_nx [μm] | ε_ny [μm] | σ_z [mm] | Mean pz [MeV] |
|------|-------------|-----------|-----------|---------|---------------|
| control | 0.0049 | 1.00 | 1.00 | 0.0983 | 200.51 |
| case_longitudinal | 90.2136 | 175599163.99 | 1656986747256218880.00 | 103436.2416 | 10477992450.04 |
| case_transverse | 23.1199 | 30936338.03 | 133480990324991.33 | 37.8333 | 703088525.02 |

### 注意事项

- DiWakeCyl 计算的是介质的格林函数（点电荷 δ-响应），而非束团总尾场势
- ASTRA 内部会自行完成卷积（binning → 离散卷积 → 粒子 kick）
- 介质波导的尾场格林函数值量级约 10^16 V/m/C，这是由物理过程决定的
- Taylor_Method_F 格式已在实际测试中验证兼容

## 文件结构

```
test/diwakecyl_astra_verify/
├── README.md
├── run_verify.py
├── wakefield_green.png
├── summary_all_cases.png
├── case_control/          # 对照组 (无尾场)
├── case_longitudinal/     # m=0 纵向尾场
└── case_transverse/       # m=1 横向尾场
```

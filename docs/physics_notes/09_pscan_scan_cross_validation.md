# 09 — PScan / Scan 文件的真实数据交叉验证

## 数据来源
本地 ASTRA 真实运行 (data/pscan_validate, data/scan_validate,
gitignored), 产物作为 golden 提交:

* examples/Cavity_Example/golden/astra.PScan.001 + astra_pscan.ref.001
  (TWS 9-cell 算例, Phase_scan=T, Auto_Phase=F, LSPCH=F)
* examples/Manual_Example/Example.Scan.001
  (&SCAN: Scan_para='MaxB(1)', S=0.1..0.5, S_numb=5,
   FOM(1)='hor. Emittance', FOM(2)='rms bunch length',
   FOM(3)='mean beam energy')

测试 (test/test_batch_c_fixes.py) 不依赖可执行文件, 只读 golden。

## PScan (手册 4.9, Table 4: phase[deg] E_kin[MeV] 压缩因子 beta/beta0)

1. **单腔余弦律**: E(φ) = E0 + A·cos(φ−φ0)。实测拟合残差
   5.6 keV, 相对振幅 10.47 MeV 的 0.05%; 峰谷差 = 2A 到 0.2%。
2. **与参考粒子对照**: PScan 相角 0 的能量 vs 同一运行 ref.001 的
   末态参考粒子动能, 吻合到 0.38%。残余差异来自扫描粒子与分布参考
   粒子之间约 2° 的注入相位偏移 (分布参考粒子带 z/clock 偏移), 属
   预期物理, 容差取 0.5%。
3. **压缩因子/速度比**: 全程有界; 峰相位处 ≈ 1.0 (无压缩)。
4. 注意: PScan 是单粒子扫描, 不含空间电荷; 束团平均能量见 Zemit。
   实测: 同一运行 Zemit 末行 E = 107.95 MeV vs PScan 峰 110.47 MeV
   (差 = 空间电荷 + 参考粒子相位偏移), 勿混淆两个量。

## Scan (手册 4.9, &SCAN namelist, Table 4: para z FOM(1..10))

与同一算例 (Manual_Example, MaxB=0.35 T) 的 golden Xemit/Zemit 对照:

* FOM(1) 'hor. Emittance' = 归一化水平发射度 [pi mm mrad]:
  扫描点 0.3/0.4 T 的 0.9999/1.0008 夹住 golden Xemit 末行 1.0003,
  在 0.35 T 处插值吻合 < 0.5%。
* FOM(2) 'rms bunch length' [mm] = 0.6595227, golden Zemit σz =
  0.65952, 逐点精确。
* FOM(3) 'mean beam energy' [MeV] = 999.989, golden Zemit E_kin =
  999.99, 逐点精确。

结论: read_scan 的列语义与 ASTRA 的 FOM 关键字一一对应, 已独立验证。

# 08 — Zemit 第7列 (correlation) 的单位与语义

## 问题
Zemit.001 第 7 列在手册 Table 4 中记为 <z·E'_avr>，单位为 **keV**
（不是 keV·mm）。早期实现按 keV·mm 处理（乘 1e3·1e-3 = 1），数值上
差 1000 倍。

## 手册依据
Astra-Manual_V3.2.md 的 Table 4（Zemit 行）:

    Zemit | z m | t ns | E_kin MeV | z_rms mm | dE_rms keV |
           eps_z,norm pi keV mm | z·E'_avr keV

注意 Xemit 的对应列是 <x·x'_avr> mrad = cov(x,x')/sigma_x
（见 notes/05），因此 Zemit 列按类比应为 cov(z, E_kin)/sigma_z，
量纲为能量（keV），不带长度。

## 数值验证（cross-validation）
对 examples/Manual_Example/Example.0150.001 的活粒子
（status > 1）直接计算:

    cov(z, E_kin) / sigma_z   [eV]   vs   文件第7列 x 1e3

二者吻合到 <0.5%。因此读取时:

    corr = data[:, 6] * 1e3      # keV -> eV

（原实现 *(1e3 * 1e-3) = *1 是错的。）

## 结论
* Zemit corr 列 = cov(z, E_kin)/sigma_z，内部 SI 单位为 eV。
* 测试: test/test_batch_a_fixes.py::test_zemit_corr_matches_particle_covariance。

# 06 Sigma 矩阵归一化之谜 (已破解)

## 已知事实
Sigma 文件 23 列: z[m], E_kin[MeV], 6x6 协方差上三角 21 元素
(手册 Table 4)。坐标系: (x, p~x, y, p~y, z, E_kin) (手册 §4.13)。

## 破解 (2026-08, 数值验证)

文件把动量列 (p~x/p~y) 归一化到 **mc**、能量列 (E_kin) 归一化到
**mc^2** (均无量纲):

    sig22_file * (mc)^2  = (sigma_x' * p_ref)^2     比值 1.000000
    sig66_file * (mc)^2  = sigma_E^2               比值 1.000000
    sig56_file * (mc)    = cov(z, E) [m.eV]        与 Zemit corr 列一致

历史"3.83 因子" = 1/(mc[MeV])^2 = 1/0.51099895^2 = 3.8297。
(与 sig22/sig66 差一个 3.83 倍只是 MeV 单位下的表象。)

## 处理
read_sigma_file 现在把矩阵统一换算到 SI (位置 m、动量 eV/c、
能量 eV), 并导出归一化 eigen-emittance (|imag(eig(Sigma J))|/mc),
与 Xemit eps_n 逐行对照 < 8% (enz 与 Zemit eps_zn < 0.01%)。
eigen-emittance 图的 "experimental" 标注已移除。

## 代码位置
astra_tools/io/astra_emit.py (read_sigma_file)
交叉验证: test/test_cross_validation.py::test_sigma_eigen_emittances

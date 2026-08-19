# 11 — IO 格式要点 (2026-08 审查确认)

手册 Table 1/3/4 之外的实测格式事实, 供读者/调试参考:

## 3D 场图头 (两种形式, 实测兼容)

* 逐值: n 后跟 n 个网格值, 自由格式可换行 (3D_test.*, 3D_Dipole.*);
* 紧凑: n, min, spacing (DESY laser.dat);
* read_3d_field_map 双分支判定: 前三行恰各 3 个 token 且剩余数据
  长度与网格积精确相等 -> 紧凑头; 否则逐值 (避免 3D_Dipole.bx
  那种每行 3 个值的换行巧合误判);
* 数据顺序 x 最快 (Fortran 序), F[ix,iy,iz]; 验证: 3D_test.bx 轴上
  Bx=0、y 反对称 100% (test_batch_a_fixes)。

## 尾场表 (多块)

* 首行 (nblocks, 0), 每块 (N, 0) + N 行 (s, W); 块 1 = 单极,
  其后为双极分量 (TESLA_MODULE_WAKE_TAYLOR.dat: 3 块, 20001 点);
* read_wake_potential 返回第一块, 全部块在 .blocks。

## track 文件 (8 列)

* seq | status | z | x | y | Ez | (Er 或 Ex) | (0 或 Ey);
  笛卡尔 3D 场下第 7/8 列是 Ex/Ey, 读器保留文件名语义 (Er/Ey)。

## 二进制分布 9/10 列歧义

* 布局: 9 列 = x y z px py pz clock charge status; 10 列 =
  x y z px py pz clock charge species status。
* 手册 Table 1: 第 9 列是**粒子种类** (1=电子, 2=正电子, 3=质子,
  4=氢离子, 5-14=用户定义), **不是**顺序递增的粒子编号 — 旧笔记
  以"第 9 列 = 顺序递增编号"为前提的消歧描述是错误前提, 已废弃。
  事实依据: 手册 Table 1 定义 + 仓库内全部真实 ASCII 分布文件核对,
  第 9 列恒为 1 (电子)。
* 消歧: 粒子数同时满足 N%9==0 与 N%10==0 (即 N 为 90 的倍数) 时,
  文件长度无法区分 9/10 列; species 列恒值, "形似顺序递增编号"的
  启发式在真实 species 列上永不命中 — 判定必须依据列语义
  (第 9 列 species 恒值 / 第 10 列 status 变化); 现行读器的消歧
  修正由并行代码任务进行。

## Sigma 归一化 (见 notes/06)

* 动量列 ÷mc、能量列 ÷mc (均无量纲; 读回 SI 时动量/能量列乘 mc,
  交叉项乘 mc 一次), 即 scale=[1, mc, 1, mc, 1, mc] — 不是
  "能量列 ÷mc^2";
* "3.83 因子" = 1/(mc[MeV])^2 = 1/0.51099895^2 仅是 MeV 单位表象;
* 读者换算 SI 后 eigen-emittance 与 Xemit 对照 <8%, enz <0.01%。

# astra-notebook 测试方案 (完整版)

> 测试总原则: **尽一切办法验证数据正确**。任何新代码必须至少通过
> 五层中的前四层; 涉及 notebook 的改动必须过第五层。所有统计量必须
> 与 ASTRA 自身输出交叉验证, 或与解析公式对照; 二者都没有的量
> (见 physics_notes/10) 必须显式标注"未验证"。

运行入口:

    .venv/bin/python -m pytest test/ -q          # 第 1-4 层
    bash test/e2e_notebooks.sh                   # 第 5 层 (需可执行文件)

---

## 第 1 层: 单元测试 (纯函数, 无外部数据)

文件: test/test_backend.py, test/test_format_input.py, test/test_batch_b_fixes.py 的一部分。

* 覆盖: namelist 读写/解析 (含多赋值、数组、括号索引)、参数元数据
  (405 参数, Nemit 的 π 单位)、Distribution.from_arrays 总电荷、
  可执行文件定位、导出 CSV/npz 格式。
* 特征: 不读 examples/ 数据、不调用可执行文件; 秒级完成。
* 新增代码要求: 每个公开纯函数至少一个用例 (含边界: 空输入、单行、
  负值、混合符号)。

## 第 2 层: 黄金样本回归 (不依赖可执行文件)

文件: test/test_golden_examples.py; 数据: examples/*/golden/ 与
examples/golden_expected.json。

* 覆盖: 9 个官方算例的归档输出 (Xemit/Zemit/ref/Log/PScan/Scan) 与
  记录在 golden_expected.json 的期望末行值比对 (rel < 0.5%)。
* 黄金样本生成流程 (仅在本地重跑后更新):
  1. notebooks/08_examples.ipynb 的 run_example() 在 data/ 下真跑;
  2. 核对新输出与 golden 的差异 (compare_xemit 打印 rel%);
  3. 确认差异为物理/版本变化后, 更新 golden 文件与
     golden_expected.json, 并在提交说明中记录原因。
* 注意: golden_expected.json 由我们自己的解析器生成, 是"自洽回归",
  独立性由第 3 层保证。

## 第 3 层: 交叉验证 vs ASTRA 自身输出 (核心层, "尽一切办法")

文件: test/test_cross_validation.py (统计链), test/test_batch_a_fixes.py
(Zemit corr / 3D 场图 / 加权 / slice), test/test_batch_c_fixes.py
(PScan / Scan)。

已建立的独立对照:

| 我们的量 | ASTRA 数据 | 容差 |
|---|---|---|
| mean/sigma_x, σx', εnx | Xemit 末行 | <0.5% (实测 <0.02%) |
| mean E_kin, σz, σE, εnz | Zemit 末行 | <0.5% |
| cov(z,E)/σz | Zemit 第 7 列 (keV×1e3) | <0.5% |
| 3D 场图轴序 | 3D_test.bx 轴上 Bx=0 + y 反对称 100% | 精确 |
| BFF | 解析高斯 exp(−(kσz)²) | rel 5% |
| slice εn | 解析高斯 σxσpx/mc² | rel 15% (每 slice) |
| PScan E(φ) | 同运行 ref.001 参考粒子 | 0.5% (相位偏移, 见 notes/09) |
| Scan FOM(1/2/3) | golden Xemit εnx / Zemit σz, E | <0.5% (实测逐列精确) |
| 核心分数曲线 f=1.0 | compute_statistics | 严格一致 |

新增统计量的要求: 优先找 ASTRA 文件对照; 没有则解析公式; 再没有则
自洽 (f=1 全束团) + 在 physics_notes/10 登记"未验证"。

## 第 4 层: 绘图正确性

文件: test/test_plots.py, test/test_batch_c_fixes.py 的绘图部分。

* 每个绘图函数断言: 轴标签含正确单位串 (如 "z [m]", "t [ns]",
  "[pi mm mrad]")、图例含物理量名、KDE 密度系统被使用。
* 黄金场图 (3D_test.*) 的真实对称性作为场图数据正确性的物理锚点。
* 新增绘图必须: (1) 标签/单位断言; (2) 若涉及新数据源, 加一层
  第 3 层的数据对照。

## 第 5 层: 端到端 (notebook 执行, 需 ASTRA/Generator 可执行文件)

文件: test/e2e_notebooks.sh。

* 逐个执行 notebooks/00..08:

      jupyter nbconvert --to notebook --execute notebooks/XX.ipynb \
        --ExecutePreprocessor.kernel_name=astra-notebook

* 08 会真跑全部 9 个官方算例并做黄金比对 (约 1-2 分钟)。
* 失败标准: 任何单元格异常或 nbconvert 返回非零。
* 本机前提: astra/generator 在 PATH 或项目 ASTRA/ 目录;
  Plasma_Example_2 需要本地 laser.dat (65MB, gitignored)。

## 测试清单 (当前 83 项)

* test_backend.py — 单元 + 导出 + 可执行文件
* test_format_input.py — namelist 生成
* test_cross_validation.py — 统计链 vs ASTRA (4)
* test_golden_examples.py — 9 算例黄金回归 (7)
* test_plots.py — 绘图标签/单位 (24)
* test_batch_a_fixes.py — 批 A 物理修复红绿测试 (8)
* test_batch_b_fixes.py — 批 B 工程修复 + BFF (12)
* test_batch_c_fixes.py — 批 C 展示 + PScan/Scan 交叉验证 (14)

## 新增功能的标准流程 (写代码前先写测试)

1. 若涉及物理量: 先在 physics_notes/ 记下 手册依据 → 数值验证计划;
2. 写红测试 (先跑, 确认失败);
3. 实现代码 → 测试转绿;
4. pyflakes 清零; 完整 pytest;
5. 若改了 notebook: 重跑 e2e;
6. 提交说明写明"验证了什么、对照了什么、容差"。

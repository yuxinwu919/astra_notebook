# astra-notebook 用户手册

astra-notebook 是 ASTRA / Generator 的前后处理工作台: 前端为 8 个
Jupyter Notebook (表单点选 + 现代绘图), 后端为 Python 包 astra_tools。
它替代 macOS 上缺失的官方 postpro / lineplot / fieldplot。

## 1. 安装

    python3 -m venv .venv && source .venv/bin/activate
    pip install -r requirements.txt
    python -m ipykernel install --prefix .venv --name astra-notebook

ASTRA / Generator 可执行文件放入系统 PATH 或项目 ASTRA/ 目录。

    jupyter notebook     # 在项目根目录启动

## 2. 标准工作流

| 步骤 | Notebook | 做什么 |
|------|----------|--------|
| 1 | 01_generator | 表单设参数 -> 运行 Generator -> 束团预览 |
| 2 | 02_astra | 环境自检 + 追踪设置 (表单或 .in 文本) + 运行 + 输出清单 |
| 3 | 03_postpro | 相空间图 + 统计表 + 切片/BFF/导出 (选 z 位置) |
| 4 | 04_lineplot | 束流参数随 z 演化 (九图 + 速度/步长 + t 轴变体) |
| 5 | 05_fieldplot | 腔场/螺线管/3D 场图/激光/等离子体 |
| 6 | 06_examples | 官方算例展示 + 一键复现与黄金比对 |

每个 Notebook 独立可运行 (复制项目文件夹到任意路径均可)。

## 3. 参数表单用法

- 表单由 ASTRA 手册第 6/7 章的参数元数据库自动生成, 悬停控件可见
  手册描述与单位;
- 02 中基础组 (NEWRUN/OUTPUT) 预填了可运行的默认值 (纯漂移), 其余
  组 (CAVITY/SOLENOID/...) 只写入你改动过的参数;
- 文本模式: 把现成 .in 复制为 data/workspace/astra.in 即可跳过表单;
- 提示: ASTRA 输入卡中的文件路径用相对路径 (工作目录), 过长的绝对
  路径会被 ASTRA 截断。

## 4. 单位与物理约定 (重要)

- 发射度: 图与表均显示 "π mm mrad", 数值与 ASTRA 打印完全一致
  (π 表示 RMS 相空间椭圆面积语义, 数值上即 mm·mrad); 纵向为
  "π keV mm" (数值 keV·mm)。
- 相空间图: x/y [mm], x'/y' [mrad], z [mm], dp/p [%]; 束流统计:
  能量 MeV/keV, 电流 A, 电荷 nC, 场 MV/m 与 T。
- 密度图默认裁剪 0.5% 极端离群点并在图上标注, 保证主体结构清晰。
- 螺线管中的发射度按手册 4.13.1 用正则动量计算 (自动处理)。

## 5. 数据导出

07 中一键导出: CSV (表头带单位注释, 适合外部绘图) 与 npz (原始数组)。

## 6. 常见问题

- 图太挤/点聚成一团: 已内置百分位裁剪; 若仍不满意可在绘图函数
  调大 clip_q 参数。
- 提示未找到可执行文件: 放入 PATH 或项目 ASTRA/ 目录。
- ASTRA 报 Program stops / 输入卡解析失败: 用 02 的表单模式重新
  生成输入卡; 检查文件路径为相对路径。
- 中文图例乱码: matplotlib 用 DejaVu 字体, 图内标签均为英文。

# astra-notebook 开发手册

## 架构

前端 (notebooks/) 只做编排与呈现; 一切数据逻辑在后端包 astra_tools/
(import astra_tools)。复制文件夹即用, 无打包。

    astra_tools/
    |-- io/        读取: 相空间(二进制9/10列+ASCII)、Xemit/Yemit/Zemit/
    |              Cemit/LandF(通用表驱动)、ref/Sigma、场图(1D/TWS/3D)
    |-- deck/      输入卡: metadata(手册提取的405参数元数据) + namelist 读写
    |-- analysis/  统计/发射度(正则动量)/切片/BFF (SI 内部单位)
    |-- plot/      KDE 密度体系 + postpro/lineplot/fieldplot 图形
    |-- run/       可执行文件定位/运行(失败标记检测)/输出发现/备份
    |-- widgets/   ipywidgets 组件 (选择器/表单/面板)
    `-- export.py  CSV/npz 导出

## 数据流

文件 -> io 读取 (SI + 绝对坐标归一化) -> analysis 计算 (SI) ->
plot 渲染 (显示单位在边界转换) -> export 导出。

## 物理铁律

见 AGENTS.md 的 Physics rules 8 条 (status>1、相对坐标、动量非动能、
正则动量、π 单位等)。任何改动须先改 test/ 中的对应验证再动代码。

## 测试 (五层)

    pytest test/                          # 1-4 层: 单元/黄金/交叉/绘图
    # 第 5 层: 逐个执行 notebooks (需可执行文件)
    jupyter nbconvert --to notebook --execute notebooks/XX.ipynb \
      --ExecutePreprocessor.kernel_name=astra-notebook

## 添加新绘图

1. 实现于 astra_tools/plot/, 输入 SI、输出显示单位;
2. 用 plot/_density.density2d 做密度, clip_percentile 做范围;
3. 必须含完整图例与单位标签;
4. 在 test/test_plots.py 增加标签/单位断言。

## 参数元数据库

astra_tools/data/parameters/*.json 由解析脚本从手册 md 提取 (脚本见
git 历史 /tmp/parse_manual_params.py); 手册更新时重跑该脚本。

## 黄金样本

examples/<name>/golden/ 由本地真跑生成 (data/golden_runs/, gitignore);
期望值在 examples/golden_expected.json; 重跑后须重新生成期望值并核对差异。

## 已知边界 (手册第 5 章的展示类末梢, 数据层均已支持)

主体功能已全部覆盖; 以下纯展示变体未实现, 数据源 (对应文件读取器)
均已就绪, 需要时按 plot/ 现有模式补充即可:

* postpro: 三视图 vs 时间 (t 轴变体)、核心束长 vs 电荷、含孔径几何
  叠加、3D 椭圆的交互旋转;
* lineplot: 粒子速度曲线 (可由 beta 推出)、平均步长曲线;
* fieldplot: 阴极表面场、激光场/等离子体场 vs z/zeta 的专用图
  (3D 场图截面已可查看相应数据)。

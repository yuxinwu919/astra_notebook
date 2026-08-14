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

## 已知边界 (第 5 章展示类末梢)

已实现:
* postpro: 三视图 vs 时间 (所有演化图的 x_axis='t' 参数)、核心束长/
  发射度 vs 电荷分数 (plot_core_fraction_curves, 自算)、含孔径几何
  叠加 (plot_envelope_with_aperture + aperture_elements);
* lineplot: 粒子速度曲线 (plot_velocity_evolution)、平均步长曲线
  (plot_step_size_evolution);
* fieldplot: 阴极表面场 (plot_cathode_emission include_spch)、激光
  3D 图轴上剖面 vs z/t (plot_laser_on_axis)、等离子体密度剖面
  (plot_plasma_profile)。

剩余 (刻意不做或留待后续):
* 3D 椭圆的交互旋转 —— 项目约定静态图优先, 已提供多角度静态视图
  (plot_slice_ellipses_3d);
* PScan/Scan/Error 已有真实数据交叉验证 (test_batch_c_fixes.py),
  但 Error 文件尚无真实 golden (无 ErrorS 算例)。

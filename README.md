# ASTRA & Generator — Jupyter 交互式模拟界面

本项目提供了一套完整的 Python 驱动工作流，用于调用 **ASTRA**（粒子追踪程序）和 **Generator**（初始分布生成程序）。

## 项目结构

项目拆分为两个独立 Notebook + 一个共享工具模块：

| 文件 | 用途 | 依赖 |
|------|------|------|
| `generator_interface.ipynb` | 初始束团分布生成 | 仅需 Generator 可执行文件 |
| `astra_interface.ipynb` | 粒子追踪模拟 & 相空间分析 | 仅需 ASTRA 可执行文件 |
| `utils.py` | 共享工具模块（解析/统计/绘图） | 两个 Notebook 共用 |
| `format_input.py` | 输入文件格式规范化工具 | 独立 CLI 工具 |

**设计优势**：
- **可独立使用**：ASTRA Notebook 可加载任意 `.ini` 分布文件，不依赖 Generator
- **可组合使用**：先用 Generator Notebook 生成分布，再用 ASTRA Notebook 追踪
- **代码复用**：所有共享逻辑集中在 `utils.py`，避免重复维护

## 设计理念

- **摒弃手工编辑输入文件**：所有模拟参数以 Python 字典形式定义，一键生成 Fortran namelist 文件
- **项目隔离**：所有输入/输出文件均存放在 Notebook 所在的项目目录，不污染 ASTRA 系统目录
- **交互式分析**：运行模拟后立即进行统计计算和数据可视化

## 环境配置

### 1. 系统要求

| 组件 | 要求 | 说明 |
|------|------|------|
| 操作系统 | macOS / Linux / Windows (WSL 推荐) | 本项目在 macOS 上开发测试 |
| Python | 3.9 或更高版本 | 推荐使用虚拟环境或 Conda |
| ASTRA | 可执行文件位于项目目录 | 粒子追踪模拟引擎 |
| Generator | 可执行文件位于项目目录 | 初始粒子分布生成器 |

### 2. Python 环境搭建

**方式一：pip + venv（推荐）**

```bash
# 1. 克隆或进入项目目录
cd astra_notebook

# 2. 创建虚拟环境
python3 -m venv venv

# 3. 激活虚拟环境
# macOS / Linux:
source venv/bin/activate
# Windows:
venv\Scripts\activate

# 4. 安装依赖
pip install -r requirements.txt

# 5. 启动 Jupyter
jupyter notebook
```

**方式二：Conda（适合科学计算用户）**

```bash
# 1. 创建环境
conda create -n astra_env python=3.10 -y

# 2. 激活环境
conda activate astra_env

# 3. 安装核心依赖
conda install numpy matplotlib scipy jupyter -y

# 4. 安装粒子解析器（二选一）
pip install ocelot
# 或
pip install lume-astra

# 5. 启动 Jupyter
jupyter notebook
```

### 3. ASTRA / Generator 可执行文件

ASTRA 和 Generator 是专有科学计算软件，可从 DESY 获取：https://www.desy.de/~mpyflo/

程序查找可执行文件的顺序为：

1. **系统 PATH** — 如已安装到 `/usr/local/bin/` 等处
2. **项目 `astra/` 子目录** — 项目自包含方案，推荐

**推荐方式（项目自包含）：**

将 `astra` 和 `generator` 可执行文件放入 `astra/` 目录：

```
astra_notebook/
└── astra/           # ASTRA 发行版目录
    ├── astra        # 当前平台可执行文件
    ├── generator
    ├── Astra_for_MacOS_AppleSilicon/  # 各平台备份
    ├── Astra_for_Windows/
    └── ...
```

> `astra/` 目录下已存放了不同平台的可执行文件。切换平台时，将对应子目录中的文件复制到 `astra/` 根目录并覆盖即可。

**macOS 安全提示：**

首次运行时，macOS Gatekeeper 可能阻止未签名程序运行。请在终端中执行：

```bash
cd astra_notebook
xattr -d com.apple.quarantine astra/astra astra/generator 2>/dev/null
chmod +x astra/astra astra/generator
```

**其他平台：**

- **Linux**：确保文件有执行权限（`chmod +x`）
- **Windows**：推荐使用 WSL2；原生运行需将 `.exe` 文件放在 `astra/` 目录

### 4. 验证环境

在项目目录下运行以下命令确认环境就绪：

```bash
python -c "
import numpy, matplotlib, scipy; 
print('✓ 核心库就绪')
"
```

启动 Jupyter 后，依次运行 `generator_interface.ipynb` 中的前两个代码单元格（「初始化」和「Generator 输入参数」）。若无报错，则环境配置完成。

### 5. 依赖说明

| 包名 | 必需 | 用途 |
|------|------|------|
| `numpy` | ✓ 必需 | 数值计算与二进制文件解析 |
| `matplotlib` | ✓ 必需 | 相空间与分布绘图 |
| `scipy` | ✓ 必需 | 统计分析（高斯拟合等） |
| `jupyter` | ✓ 必需 | Notebook 运行环境 |
| `ocelot` | 推荐 | 首选粒子分布解析器，功能最全 |
| `lume-astra` | 备选 | 轻量级纯 Python 解析器 |
| ASTRA 可执行文件 | ✓ 必需 | 粒子追踪模拟 |
| Generator 可执行文件 | ✓ 必需 | 初始分布生成 |

> **注意**：即使 `ocelot` 和 `lume-astra` 都未安装，Notebook 仍可使用内置的 `np.fromfile` 二进制解析器正常运行。但推荐至少安装一个以获得更好的数据兼容性。

## 工作流概览

### 方式一：组合使用（完整流程）

1. 打开 `generator_interface.ipynb` → 设置 `gen_params` → 运行 → 输出 `bunch.ini`
2. 打开 `astra_interface.ipynb` → `Distribution` 自动指向 `bunch.ini` → 运行 → 追踪结果

### 方式二：独立使用 ASTRA

直接打开 `astra_interface.ipynb`，将 `DISTRIBUTION_FILE` 设为任意 `.ini` 分布文件路径即可进行追踪模拟，无需运行 Generator。

## 快速开始

1. 确保已完成 [环境配置](#环境配置) 所有步骤
2. 在项目根目录启动 Jupyter Notebook：
   ```bash
   cd astra_notebook
   jupyter notebook
   ```
3. **生成初始分布**：打开 `generator_interface.ipynb` → 依次运行所有单元格 → 得到 `bunch.ini`
4. **粒子追踪模拟**：打开 `astra_interface.ipynb` → 依次运行所有单元格 → 得到追踪结果和相空间分析
5. 修改参数只需编辑对应字典（`gen_params`、`astra_newrun_params` 等），重新运行相关单元格

## 常见问题

### Q: 两个 Notebook 有什么区别？我应该用哪个？
**A:** 
- `generator_interface.ipynb`：生成初始粒子分布（`bunch.ini`），仅需 Generator
- `astra_interface.ipynb`：对已有分布文件进行粒子追踪模拟，仅需 ASTRA
- 如果你已有 `.ini` 分布文件，直接使用 `astra_interface.ipynb` 即可

### Q: 两个 Notebook 可以独立运行吗？
**A:** 可以。`astra_interface.ipynb` 的 `DISTRIBUTION_FILE` 变量可指向任意 `.ini` 文件，不依赖 Generator Notebook。

### Q: 提示 "未找到可执行文件 'astra'"
**A:** 请确保 `astra` 和 `generator` 可执行文件已正确放置。参见 [ASTRA / Generator 可执行文件](#3-astra--generator-可执行文件) 章节。

### Q: ocelot 安装失败
**A:** ocelot 依赖较复杂，尤其在 Windows 上。替代方案：
```bash
pip uninstall ocelot
pip install lume-astra
```
或两个都不安装，Notebook 将自动使用内置解析器。

### Q: 中文图表显示方块
**A:** 系统缺少中文字体。编辑 Notebook 第 1 个代码单元格中的字体设置：
```python
plt.rcParams['font.sans-serif'] = ['DejaVu Sans']  # 使用英文字体
```
或安装中文字体（macOS: 系统自带；Linux: `sudo apt install fonts-wqy-microhei`）。

### Q: Generator 输出文件未找到
**A:** Generator 可能使用了不同的默认输出文件名。Notebook 会自动在工作目录中搜索 `*.ini` 和 `*.dat` 文件。检查 `simulation_files/` 目录下的文件。

## 文件结构

```
astra_notebook/
├── generator_interface.ipynb         # Generator Notebook（初始分布生成）
├── astra_interface.ipynb             # ASTRA Notebook（粒子追踪模拟，8 个 namelist）
├── utils.py                          # 共享工具模块（解析/统计/绘图/namelist 生成）
├── format_input.py                   # 输入文件格式规范化 CLI 工具
├── README.md                         # 项目说明文档
├── requirements.txt                  # Python 依赖清单
├── backup.py                         # 仿真结果备份脚本
├── astra/                            # ASTRA 发行版（含多平台执行文件）
│   ├── astra                         # ASTRA 可执行文件
│   └── generator                     # Generator 可执行文件
├── examples/                         # ASTRA 官方示例
│   ├── Manual_Example/               # 手册示例
│   ├── Wake/                         # 尾场算例
│   └── Aperture/                     # 孔径算例
├── docs/                             # 文档
│   └── astra_wakefield_guide/        # 尾场机制详解（LaTeX）
├── test/                             # 测试数据
│   └── WHPS_LINAC_test_data/         # WHPS LINAC 实际算例
└── simulation_files/                 # 默认工作目录（自动创建）
    ├── generator.in                  # Generator 输入文件
    ├── bunch.ini                     # Generator 输出 / ASTRA 输入
    ├── astra.in                      # ASTRA 输入文件
    └── astra.*.001                   # ASTRA v4 输出文件
```

## 参数说明

### Generator (`&INPUT` namelist)

所有参数定义在 `gen_params` 字典中，完整参数列表参考 ASTRA 手册第 7 章。

| 参数 | 含义 | 单位 | 默认值 |
|------|------|------|--------|
| `FNAME` | 输出分布文件名 | — | `'bunch.ini'` |
| `IPart` | 生成的粒子数 | — | 10000 |
| `Species` | 粒子种类 | — | `'electrons'` |
| `Q_total` | 束团总电荷 | nC | 0.1 |
| `Probe` | 是否生成探针粒子 | T/F | T |
| `Noise_reduc` | Hammersley 准随机降噪 | T/F | T |
| `Cathode` | 阴极发射模式（T=时间分布, F=空间分布） | T/F | T |
| `Ref_Ekin` | 参考粒子动能 | MeV | 0.0 |
| `Ref_zpos` | 参考粒子 z 位置 | m | 0.0 |
| `Dist_z` | 纵向分布类型 | — | `'uniform'` |
| `sig_z` | RMS 束团长度 | mm | 0.0 |
| `C_sig_z` | Gauss 截断倍数 | — | 0.0 |
| `Dist_pz` | 纵向动量分布类型 | — | `'uniform'` |
| `sig_Ekin` | RMS 能量散度 | keV | 0.0 |
| `cor_Ekin` | z–E 相关能散 | keV | 0.0 |
| `Dist_x`, `Dist_y` | 横向分布类型 | — | `'gaussian'` |
| `sig_x`, `sig_y` | RMS 横向尺寸 | mm | 1.0 |
| `Dist_px`, `Dist_py` | 横向动量分布类型 | — | `'gaussian'` |
| `Nemit_x`, `Nemit_y` | 归一化横向发射度 | π·mrad·mm | 0.0 |
| `cor_px`, `cor_py` | 相关散角 | mrad | 0.0 |

**注意**：
- `Cathode=F` 时必须设置 `sig_z`（空间分布）；`Cathode=T` 时必须设置 `sig_clock`（时间分布）
- 字符串值需包含单引号（如 `"'gauss'"`）
- 分布式关键字仅首字符有效：`'g'`=`'gauss'`, `'u'`=`'uniform'`, `'p'`=`'plateau'`

### ASTRA namelists（所有手册参数，按需开关）

ASTRA Notebook 支持 **8 个 namelist 模块**，参数均按手册顺序排列，通过开关控制是否写入 `astra.in`：

| Namelist | 手册章节 | 开关变量 | 默认 | 参数数 |
|----------|----------|----------|------|--------|
| `&NEWRUN` | §6.1 | 始终启用 | — | 13 |
| `&OUTPUT` | §6.2 | 始终启用 | — | 35 |
| `&CHARGE` | §6.6 | `USE_CHARGE` | `True` | 28 |
| `&SCAN` | §6.3 | `USE_SCAN` | `False` | 18 |
| `&CAVITY` | §6.9 | `USE_CAVITY` | `False` | 28 |
| `&SOLENOID` | §6.10 | `USE_SOLENOID` | `False` | 12 |
| `&WAKE` | §6.8 | `USE_WAKE` | `False` | 17 |
| `&APERTURE` | §6.7 | `USE_APERTURE` | `False` | 9 |

**两种输入模式**：
- **模式 A**（`USE_EXISTING_ASTRA_INPUT=True`）：导入现有 `.in` 文件，直接预览并运行
- **模式 B**（`USE_EXISTING_ASTRA_INPUT=False`）：通过 Python 字典自动生成 `astra.in`

**重要**：`Distribution` 必须与 Generator 的输出文件名一致。

## 输入文件格式规范化工具

`format_input.py` 用于将 ASTRA / Generator 输入文件统一为规范格式：

```bash
# 原地格式化单个文件
python format_input.py examples/Aperture/astra.in

# 仅检查不修改
python format_input.py examples/Wake/Wake_Files/Wake.in --check

# 输出到指定文件
python format_input.py examples/Manual_Example/Example.in -o simulation_files/astra.in

# 批量格式化项目中所有 .in 文件
python format_input.py --all
```

**格式化规则**：
| 规则 | 说明 |
|------|------|
| CRLF → LF | Windows 换行符自动转换 |
| 行尾空白 | 去除多余空格和 Tab |
| key=value | 统一等号前后无空格，行尾加逗号 |
| 缩进 | 统一 2 空格 |
| 块间空行 | 多余空行压缩为一空行 |
| `&` / ` /` | 统一块开始/结束标记格式 |

> **提示**：Windows 上创建的 `.in` 文件通常含 CRLF 换行符，macOS/Linux 上的 Fortran 运行时无法正确解析。在导入前建议运行 `python format_input.py --all` 批量处理。

## 数据解析器优先级

1. **ocelot**（首选）：`astraBeam2particleArray` → `ParticleArray`，支持统计计算
2. **lume-astra**（备选）：`lume_astra.read_astra()`
3. **自定义 numpy 解析器**（回退）：`np.fromfile` 直接读取二进制格式

## 二进制文件格式（ASTRA 分布文件）

参考 ASTRA 手册 Table 1：

- **文件头**：5 个 float64（参考时间 [ns], 参考能量 [eV/c], 总电荷 [nC], 保留, 保留）
- **每粒子**：9 个 float64（x [m], y [m], z [m], px [eV/c], py [eV/c], pz [eV/c], clock [ns], macro_charge [nC], status_flag）
- ASCII 格式额外包含 particle_index 列（I4 整数）

## 统计量说明

- **几何发射度**：$`\varepsilon_{\text{geom}} = \sqrt{\langle x^2\rangle\langle x'^2\rangle - \langle xx'\rangle^2}`$
- **归一化发射度**：$`\varepsilon_n = \beta\gamma \cdot \varepsilon_{\text{geom}}`$
- **Twiss 参数**：$`\beta = \langle x^2\rangle/\varepsilon`$, $`\alpha = -\langle xx'\rangle/\varepsilon`$

## 可扩展方向

- **启用 CAVITY / WAKE / SCAN**：将对应 `USE_*` 开关设为 `True`，填入参数即可
- **数组参数两种写法**：
  - 列表：`"File_Efield": ["'a.dat'", "'b.dat'"]` → `File_Efield(1)='a.dat', File_Efield(2)='b.dat'`
  - 带索引键：`"File_Efield(1)": "'a.dat'", "Nue(1)": 2.9985` → 同格式输出
- **批量参数扫描**：使用 `for` 循环或 `itertools.product` 遍历多组参数
- **并行计算**：通过 `mpirun -np N astra` 使用 MPI 加速
- **自动优化**：结合 `scipy.optimize` 进行束团匹配
- **结果数据库**：使用 pandas 或 SQLite 存储多次运行的结果

## 参考

- ASTRA 用户手册 V3.2（DESY, 2017）
- Ocelot 文档：https://github.com/ocelot-collab/ocelot
- Lume-astra：https://github.com/radiasoft/lume-astra

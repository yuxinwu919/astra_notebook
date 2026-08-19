# examples/Aperture — 孔径: 圆孔 + 圆柱堵块

二轮审计 (R2-3-1) 重建后的算例布局与运行说明。

## 物理布局 (重建版, 2026-08)

输入: `astra.in` + `test.ini`(500 粒子: 1 参考粒子 + 6 探针 + 493 核心)。

| 元件 | 位置 | 几何 | 作用 |
| ---- | ---- | ---- | ---- |
| 圆孔 | 0.10–0.12 m | R=1.5 mm, 水平偏移 +0.5 mm | 刮削满足 sqrt((x−0.5mm)²+y²) > 1.5 mm 的粒子 |
| 堵块 | 0.19–0.20 m | 负半径 (圆柱 plug) | 拦停进入该区间的全部粒子 |
| ZSTOP | 0.17 m | — | 默认运行止于圆孔之后、堵块之前 |

`test.ini` 探针 (flag 3, 轨迹探针, 相对参考粒子的 mm 级偏移):

| 探针 | x [mm] | y [mm] | 圆孔处 r_eff [mm] | 结局 |
| ---- | ------ | ------ | ----------------- | ---- |
| P1 | +2.2 | 0.0 | 1.70 | 被圆孔刮削 (flag -17) |
| P2 | −2.0 | 0.0 | 2.50 | 被圆孔刮削 (flag -17) |
| P3 | 0.0 | +2.0 | 2.06 | 被圆孔刮削 (flag -17) |
| P4 | +0.5 | 0.0 | 0.00 | 穿圆孔 → 到达 ZSTOP |
| P5 | −0.4 | +0.3 | 0.95 | 穿圆孔 → 到达 ZSTOP |
| P6 | 0.0 | −0.5 | 0.71 | 穿圆孔 → 到达 ZSTOP |

核心 (493 粒子) |x|,|y| ≤ 0.3 mm → 全部穿圆孔。参考粒子在 z=0 (束团中心),
pz=250.51 MeV/c (250 MeV, 与原算例一致), 不再落在 ZSTOP 之后。

## 运行

    # 方式 1: 共享规格 (推荐, 与 notebook 一致)
    .venv/bin/python - <<'PY'
    import sys; sys.path.insert(0, ".")
    from examples._examples_spec import run_example, compare_xemit
    work = run_example("Aperture")     # 舞台化 + 跑 astra
    compare_xemit("Aperture", work)    # 护栏 + 末行比对
    PY

    # 方式 2: 手动 (需 astra 在 PATH 或项目 ASTRA/)
    cd examples/Aperture
    astra astra.in

预期 stdout: 998 迭代步、"particles lost on aperture = 3"、
"final checkpoint at z = 0.1700 m"、Xemit 500 行且末行 z=0.17 m。
预期相空间末 dump `astra.0017.001`: 500 总 / 497 active / 3 lost (flag -17,
trajectory probe lost on aperture)。sigma_x 在圆孔处从 ~158 μm 降到 ~86 μm。

## 堵块演示 (可选)

把 `astra.in` 的 `ZSTOP` 改回 0.25 m 重跑: 穿圆孔的 497 粒子在 0.19–0.20 m
被堵块全部拦停, 参考粒子止于 0.1901 m, stdout 报 "particles lost on aperture
= 500", Xemit 停在 0.19 m 之前 (末行 z 不再到达 ZSTOP — 属预期, 此配置不用于
golden 比对)。

## golden

`golden/` 由本机真实 ASTRA 二进制 (macOS Apple Silicon 构建) 确定性重跑生成
(两次运行字节一致), 包括 Xemit/Yemit/Zemit/Sigma/Log/ref 与相空间 dump
`astra.0017.001`。期望值在 `examples/golden_expected.json` 与
`test/test_golden_examples.py` (Aperture 条目, 二轮审计加入)。

## 运行卫生 (R2-3-4): IEEE_OVERFLOW_FLAG

macOS Apple Silicon 构建每次运行结束打印
"Note: The following floating-point exceptions are signalling: IEEE_OVERFLOW_FLAG"。
这是该构建的运行时噪声 (浮点异常标志未被清零), 不影响输出字节复现:
两次运行全部输出文件 (Xemit/Yemit/Zemit/Sigma/Log/ref/相空间 dump) 逐字节一致。
无需处理。

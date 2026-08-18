#!/usr/bin/env bash
# 第 5 层端到端测试: 逐个执行全部 Notebook (需 ASTRA/Generator 可执行文件)
# 覆盖 5 个任务式 notebook + 8 个单算例示例 notebook + 6 个功能演示 demo。
# 用法: bash test/e2e_notebooks.sh
set -u
cd "$(dirname "$0")/.."
VENV="$PWD/.venv"
KERNEL=${KERNEL:-astra-notebook}
PYTHON=${PYTHON:-"$VENV/bin/python"}
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/mplcfg}
# 隔离用户 shell 环境 (如激活的 anaconda/conda): 内核只从项目 venv
# 查找, PYTHONPATH 清空, 避免混入其他 Python 发行版的 site-packages
export JUPYTER_PATH="$VENV/share/jupyter"
unset PYTHONPATH

NOTEBOOKS="
notebooks/generator.ipynb
notebooks/astra.ipynb
notebooks/postpro.ipynb
notebooks/lineplot.ipynb
notebooks/fieldplot.ipynb
examples/Manual_Example.ipynb
examples/Aperture.ipynb
examples/Wake.ipynb
examples/Cavity_Example.ipynb
examples/Curved_Cathode_Example.ipynb
examples/90deg_bend_Example.ipynb
examples/Plasma_Example_1.ipynb
examples/Plasma_Example_2.ipynb
examples/postpro_demo.ipynb
examples/generator_demo.ipynb
examples/bff_demo.ipynb
examples/stats_validation_demo.ipynb
examples/lineplot_demo.ipynb
examples/fieldplot_demo.ipynb
"

FAIL=0
for nb in $NOTEBOOKS; do
  echo "=== $nb ==="
  $PYTHON -m jupyter nbconvert --to notebook --execute "$nb" \
    --output "/tmp/e2e_$(basename "$nb" .ipynb).ipynb" \
    --ExecutePreprocessor.kernel_name="$KERNEL" \
    --ExecutePreprocessor.timeout=600 >/dev/null 2>&1
  rc=$?
  if [ $rc -eq 0 ]; then
    echo "    PASS"
  else
    echo "    FAIL (rc=$rc)"
    FAIL=1
  fi
done

exit $FAIL

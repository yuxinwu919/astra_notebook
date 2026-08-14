#!/usr/bin/env bash
# 第 5 层端到端测试: 逐个执行全部 Notebook (需 ASTRA/Generator 可执行文件)
# 用法: bash test/e2e_notebooks.sh
set -u
cd "$(dirname "$0")/.."
KERNEL=${KERNEL:-astra-notebook}
PYTHON=${PYTHON:-.venv/bin/python}
export MPLCONFIGDIR=${MPLCONFIGDIR:-/tmp/mplcfg}

FAIL=0
for nb in 00_workspace 01_generator 02_astra_setup 03_run 04_postpro \
          05_lineplot 06_fieldplot 07_analysis 08_examples; do
  echo "=== $nb ==="
  $PYTHON -m jupyter nbconvert --to notebook --execute "notebooks/$nb.ipynb" \
    --output "/tmp/e2e_$nb.ipynb" \
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

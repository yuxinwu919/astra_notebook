#!/usr/bin/env bash
# 第 5 层端到端测试: 逐个执行全部 Notebook (需 ASTRA/Generator 可执行文件)
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

FAIL=0
for nb in 01_generator 02_astra 03_postpro \
          04_lineplot 05_fieldplot 06_examples; do
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

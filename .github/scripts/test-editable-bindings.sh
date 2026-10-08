#!/usr/bin/env bash
set -euo pipefail
python -m pip install build 'nanobind==2.15.0' 'scikit-build-core>=1.0' numpy scipy matplotlib lxml h5py
if [[ "$RUNNER_OS" == Linux ]]; then python -m pip install patchelf; fi
python -m pip install --no-build-isolation --no-deps -e python/pyalps --config-settings "build-dir=$PWD/_build/python"
PYALPS_TEST_DOWNSTREAM_EXPORT=1 python .github/scripts/run_with_timeout.py 900 \
  python -m pytest tests/pyalps tests/cmake -v --capture=tee-sys --junitxml=_build/ci/python.xml

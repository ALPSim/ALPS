#!/usr/bin/env bash
set -euo pipefail
python -m pip install build 'nanobind==2.15.0' 'scikit-build-core>=1.0' numpy scipy matplotlib lxml h5py
if [[ "$RUNNER_OS" == Linux ]]; then python -m pip install patchelf; fi
python -m pip install --no-build-isolation --no-deps -e python/pyalps --config-settings "build-dir=$PWD/_build/python"
suites=(python/pyalps/tests tests/cmake)
runner=(python -m pytest)
if [[ "${ALPS_CI_MPI:-0}" == 1 ]]; then
  # Initialize MPI before native bindings load the MPI-enabled runtime.
  runner=(python -c 'from mpi4py import MPI; import pytest; raise SystemExit(pytest.main())')
fi
PYALPS_TEST_DOWNSTREAM_EXPORT=1 python .github/scripts/run_with_timeout.py 900 \
  "${runner[@]}" "${suites[@]}" -v --capture=tee-sys --junitxml=_build/ci/python.xml

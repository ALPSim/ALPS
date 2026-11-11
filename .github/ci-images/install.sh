#!/bin/bash
# Installs the CI build environment inside a ci-images Dockerfile build (as root).
# PACKAGES: extra distribution packages (the matrix compilers); BOOST: a pinned Boost
# version to build into /opt/boost when the distribution's is too old.
set -euxo pipefail
. /etc/os-release
python=python3
venv_args=()
pip_packages=(numpy h5py lxml packaging 'PyYAML>=6,<7' 'jsonschema>=4,<5' pytest)
case $ID in
  ubuntu|debian)
    export DEBIAN_FRONTEND=noninteractive
    boost_pkg=libboost-all-dev
    [[ -z $BOOST ]] || boost_pkg=
    apt-get update
    # shellcheck disable=SC2086
    apt-get install -y --no-install-recommends \
      ca-certificates curl git tar zstd xz-utils file build-essential gfortran ccache pkg-config \
      python3 python3-venv python3-dev libhdf5-dev libblas-dev liblapack-dev libsqlite3-dev \
      libeigen3-dev libxml2-utils xsltproc $boost_pkg $PACKAGES
    rm -rf /var/lib/apt/lists/*
    ;;
  rocky|almalinux|rhel|centos)
    dnf install -y dnf-plugins-core epel-release
    dnf config-manager --set-enabled crb
    boost_pkg=boost-devel
    [[ -z $BOOST ]] || boost_pkg=
    python=python3.12
    # shellcheck disable=SC2086
    dnf install -y --setopt=install_weak_deps=False \
      ca-certificates git tar zstd xz file gcc gcc-c++ gcc-gfortran make ccache pkgconf \
      python3.12 python3.12-devel hdf5-devel blas-devel lapack-devel sqlite-devel \
      eigen3-devel libxslt libxml2 $boost_pkg $PACKAGES
    dnf clean all
    ;;
  fedora)
    boost_pkg=boost-devel
    [[ -z $BOOST ]] || boost_pkg=
    # Fedora's Python is often newer than PyPI's binary wheels (h5py): take the Python
    # packages from Fedora and only the build tools from PyPI.
    # shellcheck disable=SC2086
    dnf install -y --setopt=install_weak_deps=False \
      ca-certificates git tar zstd xz file gcc gcc-c++ gcc-gfortran make ccache pkgconf \
      python3 python3-devel hdf5-devel flexiblas-devel sqlite-devel \
      eigen3-devel libxslt libxml2 python3-numpy python3-h5py python3-lxml python3-packaging \
      python3-pyyaml python3-jsonschema python3-pytest $boost_pkg $PACKAGES
    dnf clean all
    venv_args=(--system-site-packages)
    pip_packages=()
    ;;
  *)
    echo "Unsupported distribution: $ID" >&2
    exit 1
    ;;
esac

"$python" -m venv "${venv_args[@]}" /opt/ci-venv
/opt/ci-venv/bin/python -m pip install --no-cache-dir --upgrade pip
/opt/ci-venv/bin/python -m pip install --no-cache-dir --only-binary=:all: \
  'cmake>=3.27' ninja "${pip_packages[@]}"
/opt/ci-venv/bin/python -c 'import numpy, h5py, lxml, yaml, jsonschema, pytest'

if [[ -n $BOOST ]]; then
  # Built with the distribution's default compiler: libstdc++'s ABI serves gcc and clang.
  PATH=/opt/ci-venv/bin:$PATH /opt/ci-venv/bin/python /tmp/ci/.github/scripts/prepare_dependencies.py "boost-$BOOST" /opt/boost
fi
mkdir -p /opt/boost

#!/usr/bin/env bash
# Build pinned Boost releases not provided by the CI runner package manager.
# Usage: CXX=g++-14 bash .github/scripts/build-boost.sh BOOST_SOURCE INSTALL_PREFIX
set -euo pipefail
source_dir="$(cd "$1" && pwd)"
mkdir -p "$2"
prefix="$(cd "$2" && pwd)"
compiler="${CXX:-c++}"
case "$compiler" in
  *clang*) toolset=clang ;;
  *) toolset=gcc ;;
esac
cd "$source_dir"
libraries=filesystem,serialization,program_options,regex,thread,date_time,chrono,test
if [[ "${ALPS_BOOST_MPI:-ON}" == ON ]]; then libraries+=,mpi; fi
./bootstrap.sh --with-toolset="$toolset" \
  --with-libraries="$libraries"
printf 'using %s : alps : "%s" ;\n' "$toolset" "$compiler" > alps-user-config.jam
if [[ "${ALPS_BOOST_MPI:-ON}" == ON ]]; then printf 'using mpi ;\n' >> alps-user-config.jam; fi
./b2 --user-config=alps-user-config.jam "toolset=$toolset-alps" \
  variant=release link=shared runtime-link=shared threading=multi cxxstd=17 \
  "-j${CMAKE_BUILD_PARALLEL_LEVEL:-2}" "--prefix=$prefix" install

# Boost's CMake configs are relocatable. Give macOS libraries equally relocatable
# install names so consumers can extract the archive into their own workspace.
if [[ "$(uname -s)" == Darwin ]]; then
  for library in "$prefix"/lib/*.dylib; do
    [[ -L "$library" ]] && continue
    install_name_tool -id "@rpath/$(basename "$library")" "$library"
    while read -r dependency; do
      if [[ "$dependency" == "$prefix"/lib/* ]]; then
        install_name_tool -change "$dependency" "@rpath/$(basename "$dependency")" "$library"
      fi
    done < <(otool -L "$library" | awk 'NR > 1 { print $1 }')
  done
fi

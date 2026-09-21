"""Locate pyalps' CMake package for downstream native extensions."""
import argparse

from . import get_cmake_dir

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--cmake-dir", action="store_true", help="print the CMake package directory")
arguments = parser.parse_args()
if arguments.cmake_dir:
    print(get_cmake_dir())
else:
    parser.print_help()

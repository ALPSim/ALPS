#  Copyright Matthias Troyer 2010.
#   Permission is hereby granted, free of charge, to any person obtaining
#   a copy of this software and associated documentation files (the “Software”),
#   to deal in the Software without restriction, including without limitation
#   the rights to use, copy, modify, merge, publish, distribute, sublicense,
#   and/or sell copies of the Software, and to permit persons to whom the
#   Software is furnished to do so, subject to the following conditions:
#  
#   The above copyright notice and this permission notice shall be included
#   in all copies or substantial portions of the Software.
#  
#   THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS
#   OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
#   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
#   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
#   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
#   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
#   DEALINGS IN THE SOFTWARE.

if (ALPS_PATCH_VISTRAILS AND NOT VISTRAILS_GIT_DIR)
    set(SEARCH_DIRS ${PROJECT_SOURCE_DIR}/../../git/vistrails ${PROJECT_SOURCE_DIR}/../vistrails )
    find_path(VISTRAILS_GIT_DIR scripts ${SEARCH_DIRS} NO_DEFAULT_PATH)
    mark_as_advanced(VISTRAILS_GIT_DIR)
endif (ALPS_PATCH_VISTRAILS AND NOT VISTRAILS_GIT_DIR)

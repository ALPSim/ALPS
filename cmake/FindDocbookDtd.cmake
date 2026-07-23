#  Copyright Synge Todo and Matthias Troyer 2009 - 2010.
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

# This module looks for fop and will define 
# DOCBOOK_DTD_FOUND and DOCBOOK_DTD_DIR 

# For fedora boxes
file(GLOB P1 FOLLOW_SYMLINKS /usr/share/sgml/docbook/xml-dtd-4.2* )

FIND_PATH(DOCBOOK_DTD_DIR
  NAMES docbookx.dtd
  PATHS /opt/local/share/xml/docbook/ ${Boost_ROOT_DIR}/tools/boostbook /usr/share/xml/docbook/schema/dtd/4.2/
  ${P1}
  PATH_SUFFIXES 4.2 docbook-dtd-4.2
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(DocbookDtd DEFAULT_MSG DOCBOOK_DTD_DIR)

MARK_AS_ADVANCED( DOCBOOK_DTD_DIR )


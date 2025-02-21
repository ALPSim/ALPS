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
# DOCBOOK_XSL_FOUND and DOCBOOK_XSL_DIR

FIND_PATH(DOCBOOK_XSL_DIR
  NAMES fo/docbook.xsl
  PATHS /opt/local/share/xsl/docbook-xsl/ /usr/share/xml/docbook/stylesheet/docbook-xsl/  /usr/share/sgml/docbook/xsl-stylesheets/
  PATH_SUFFIXES ${subdirs}
)

#MESSAGE("BOST ROOT DIR ${Boost_ROOT_DIR}")

if (NOT DOCBOOK_XSL_DIR)
    file(GLOB subdirs RELATIVE ${Boost_ROOT_DIR}/tools/boostbook ${Boost_ROOT_DIR}/tools/boostbook/docbook-xs*)
    FIND_PATH(DOCBOOK_XSL_DIR
      NAMES fo/docbook.xsl
      PATHS ${Boost_ROOT_DIR}/tools/boostbook
      PATH_SUFFIXES ${subdirs}
    )
endif (NOT DOCBOOK_XSL_DIR)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(DocbookXsl DEFAULT_MSG DOCBOOK_XSL_DIR)

MARK_AS_ADVANCED( DOCBOOK_XSL_DIR )


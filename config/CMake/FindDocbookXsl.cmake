#  Copyright Synge Todo and Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# This module looks for fop and will define 
# DOCBOOK_XSL_FOUND and DOCBOOK_XSL_DIR

FIND_PATH(DOCBOOK_XSL_DIR
  NAMES fo/docbook.xsl
  PATHS ${Boost_ROOT_DIR}/tools/boostbook
  PATH_SUFFIXES docbook-xsl-1.73.2
)

INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(DocbookXsl DEFAULT_MSG DOCBOOK_XSL_DIR)

MARK_AS_ADVANCED( DOCBOOK_XSL_DIR )


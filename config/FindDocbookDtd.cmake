#  Copyright Synge Todo and Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

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


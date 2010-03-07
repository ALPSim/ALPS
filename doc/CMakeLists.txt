#  Copyright Synge Todo and Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

find_package(Xsltproc REQUIRED)
find_package(Doxygen REQUIRED)
find_package(Java REQUIRED)

set (JAVA_OPTS "-Xmx128000000" CACHE STRING "Runtime options for Java")
set (FOP_OPTS "-Xmx128000000" CACHE STRING "Runtime options for FOP")
set(BOOST_DIR ${Boost_ROOT_DIR})
set(top_srcdir ${PROJECT_SOURCE_DIR})
set(docdir ${PROJECT_BINARY_DIR}/doc)
set(SRC_DIR ${PROJECT_SOURCE_DIR})
set(INCLUDE_PATH ${Boost_ROOT_DIR})
set(JAVA ${JAVA_RUNTIME})

configure_file(autodoc.in ${PROJECT_BINARY_DIR}/doc/autodoc)
configure_file(bb2db.in ${PROJECT_BINARY_DIR}/doc/bb2db)
#configure_file(catalog.xml.in ${PROJECT_BINARY_DIR}/doc/catalog.xml)
configure_file(db2html.in ${PROJECT_BINARY_DIR}/doc/db2html)
configure_file(db2pdf.in ${PROJECT_BINARY_DIR}/doc/db2pdf)
configure_file(doxygen.conf.in ${PROJECT_BINARY_DIR}/doc/doxygen.conf)
configure_file(singledoc.in ${PROJECT_BINARY_DIR}/doc/singledoc)

MAKE_DIRECTORY( ${PROJECT_BINARY_DIR}/doc/alps/doxygen)

install(DIRECTORY doc COMPONENT documentation DESTINATION doc
        FILES_MATCHING PATTERN "*.html" PATTERN "*.png"
        PATTERN ".svn" EXCLUDE PATTERN "CVS" EXCLUDE PATTERN "*.xml" EXCLUDE 
        PATTERN "*.xsl" EXCLUDE PATTERN "*.qbk" EXCLUDE)


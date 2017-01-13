#  Copyright Synge Todo and Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if (NOT WIN32)
    find_package(Xsltproc REQUIRED)
    find_package(Doxygen REQUIRED)
    find_package(Java REQUIRED)
    find_package(DocbookDtd REQUIRED)
    find_package(DocbookXsl REQUIRED) 
    find_package(Fop)

    set (JAVA_OPTS "-Xmx256000000" CACHE STRING "Runtime options for Java")
    set (FOP_OPTS "-Xmx256000000" CACHE STRING "Runtime options for FOP")
    mark_as_advanced(JAVA_OPTS)
    mark_as_advanced(FOP_OPTS)

    set(BOOST_DIR ${Boost_ROOT_DIR})
    set(top_srcdir ${PROJECT_SOURCE_DIR})
    set(docdir ${PROJECT_BINARY_DIR}/doc)
    set(SRC_DIR ${PROJECT_SOURCE_DIR})
    set(INCLUDE_PATH ${Boost_ROOT_DIR})
    set(JAVA ${JAVA_RUNTIME})
    

    configure_file(autodoc.in ${PROJECT_BINARY_DIR}/doc/autodoc)
    configure_file(bb2db.in ${PROJECT_BINARY_DIR}/doc/bb2db)
    configure_file(catalog.xml.in ${PROJECT_BINARY_DIR}/doc/catalog.xml)
    configure_file(db2html.in ${PROJECT_BINARY_DIR}/doc/db2html)
    configure_file(db2pdf.in ${PROJECT_BINARY_DIR}/doc/db2pdf)
    configure_file(doxygen.conf.in ${PROJECT_BINARY_DIR}/doc/doxygen.conf)
    configure_file(singledoc.in ${PROJECT_BINARY_DIR}/doc/singledoc)
    configure_file(src/ietl.qbk.in ${SRC_DIR}/doc/src/ietl.qbk)
    
    make_directory(${PROJECT_BINARY_DIR}/doc/src)
    make_directory(${PROJECT_BINARY_DIR}/doc/html)
    make_directory(${PROJECT_BINARY_DIR}/doc/html/search)
    make_directory(${PROJECT_BINARY_DIR}/doc/html/images)
    make_directory( ${PROJECT_BINARY_DIR}/doc/alps/doxygen)

    set(ALPS_AUTODOC_SOURCES alea lattice model osiris parser scheduler random)

    set(ALPS_QBK_SOURCES hdf5 alea general ietl lattice
        model osiris overview parser scheduler tools random)
        
    set(ALPS_QBK_DEPENDS src/alps.qbk )

    file(COPY ${PROJECT_SOURCE_DIR}/doc/src DESTINATION .)
    foreach (S ${ALPS_QBK_SOURCES})
        add_custom_command(OUTPUT src/${S}.xml 
                         COMMAND echo ${S}
                         COMMAND cmake -E copy ${PROJECT_SOURCE_DIR}/doc/src/${S}.qbk src
                         COMMAND ${PROJECT_BINARY_DIR}/tool/quickbook --input-file src/${S}.qbk --output-file src/${S}.xml
                         DEPENDS quickbook src/${S}.qbk)
      set(ALPS_QBK_DEPENDS ${ALPS_QBK_DEPENDS} src/${S}.xml)
    endforeach (S)

    file(GLOB HEADERS ${PROJECT_SOURCE_DIR}/src/alps/*.hpp ${PROJECT_SOURCE_DIR}/src/alps/*.h)
    add_custom_command(OUTPUT src/general.autodoc 
                       COMMAND ./autodoc src/general.autodoc ${HEADERS}
                       DEPENDS ${HEADERS} autodoc.in)
    set(ALPS_QBK_DEPENDS ${ALPS_QBK_DEPENDS} src/general.autodoc)

    set(HEADERS ${PROJECT_SOURCE_DIR}/src/alps/hdf5.hpp)
    add_custom_command(OUTPUT src/hdf5.autodoc 
                       COMMAND ./autodoc src/hdf5.autodoc ${HEADERS}
                       DEPENDS ${HEADERS} autodoc.in)
    set(ALPS_QBK_DEPENDS ${ALPS_QBK_DEPENDS} src/hdf5.autodoc)

    file(GLOB HEADERS ${PROJECT_SOURCE_DIR}/src/ietl/*.h ${PROJECT_SOURCE_DIR}/src/ietl/*.hpp ${PROJECT_SOURCE_DIR}/src/ietl/interface/*.h ${PROJECT_SOURCE_DIR}/src/ietl/interface/*.hpp)
    add_custom_command(OUTPUT src/ietl.autodoc 
                       COMMAND ./autodoc src/ietl.autodoc ${HEADERS}
                       DEPENDS ${HEADERS} autodoc.in)
    set(ALPS_QBK_DEPENDS ${ALPS_QBK_DEPENDS} src/ietl.autodoc)

    foreach (S ${ALPS_AUTODOC_SOURCES})
      file(GLOB_RECURSE HEADERS ${PROJECT_SOURCE_DIR}/src/alps/${S}/*.hpp ${PROJECT_SOURCE_DIR}/src/alps/${S}/*.h)
      add_custom_command(OUTPUT src/${S}.autodoc 
                         COMMAND ./autodoc src/${S}.autodoc ${HEADERS}
                         DEPENDS ${HEADERS} autodoc.in)
      set(ALPS_QBK_DEPENDS ${ALPS_QBK_DEPENDS} src/${S}.autodoc)
    endforeach (S)

    add_custom_command(OUTPUT src/alps.xml 
                       COMMAND cmake -E copy ${PROJECT_SOURCE_DIR}/doc/src/alps.qbk src
                       COMMAND ${PROJECT_BINARY_DIR}/tool/quickbook --input-file src/alps.qbk --output-file src/alps.xml
                       DEPENDS quickbook ${ALPS_QBK_DEPENDS})

    add_custom_command(OUTPUT src/alps.docbook 
                       COMMAND ./bb2db src/alps.xml src/alps.docbook
                       DEPENDS src/alps.xml ${ALPS_QBK_DEPENDS} bb2db.in)

    file(COPY src/images DESTINATION html)
    if (EXISTS ${Boost_ROOT_DIR}/doc/html/boostbook.css)
      set (CSSDIR ${Boost_ROOT_DIR}/doc/html/)
    endif(EXISTS ${Boost_ROOT_DIR}/doc/html/boostbook.css)
    if (EXISTS ${Boost_ROOT_DIR}/doc/src/boostbook.css)
      set (CSSDIR ${Boost_ROOT_DIR}/doc/src/)
    endif(EXISTS ${Boost_ROOT_DIR}/doc/src/boostbook.css)

    add_custom_command(OUTPUT html/index.html 
                       COMMAND cmake -E copy ${PROJECT_SOURCE_DIR}/doc/src/alps_logo.gif html
                       COMMAND cmake -E copy ${CSSDIR}/boostbook.css html
                       COMMAND ./db2html src/alps.docbook html/
                       DEPENDS 
                       src/alps.docbook 
                       ${CSSDIR}/boostbook.css
                       ${PROJECT_SOURCE_DIR}/doc/src/alps_logo.gif
                       db2html.in
                       )


    add_custom_command(OUTPUT alps.pdf 
                       COMMAND ./db2pdf src/alps.docbook alps.pdf
                       DEPENDS 
                       src/alps.docbook 
                       db2pdf.in
                       )

    #if(FOP_FOUND)
    #  add_custom_target(docs ALL DEPENDS html/index.html alps.pdf)
    #else(FOP_FOUND)
      add_custom_target(docs ALL DEPENDS html/index.html)
    #endif(FOP_FOUND)
endif (NOT WIN32)

if (ALPS_BUILD_SOURCE)

install(DIRECTORY . DESTINATION alps/doc
        FILES_MATCHING 
        PATTERN "*.html" 
        PATTERN "*.pdf" 
        PATTERN "*.png"
        PATTERN "*.pdf" 
        PATTERN ".svn" EXCLUDE 
        PATTERN "CVS" EXCLUDE
        )

install(DIRECTORY ${PROJECT_BINARY_DIR}/doc DESTINATION alps
        FILES_MATCHING 
        PATTERN "*.pdf" 
        PATTERN "*.html" 
        PATTERN "*.png" 
        PATTERN "*.gif" 
        PATTERN "*.xsl" 
        PATTERN "*.css"
        PATTERN "CMakeFiles" EXCLUDE
        PATTERN "src" EXCLUDE
        )

else (ALPS_BUILD_SOURCE)
install(DIRECTORY doc COMPONENT documentation DESTINATION .
        FILES_MATCHING 
        PATTERN "*.html" 
        PATTERN "*.pdf" 
        PATTERN "*.png"
        PATTERN "*.pdf" 
        PATTERN ".svn" EXCLUDE 
        PATTERN "CVS" EXCLUDE
        )

install(DIRECTORY ${PROJECT_BINARY_DIR}/doc COMPONENT documentation DESTINATION .
        FILES_MATCHING 
        PATTERN "*.pdf" 
        PATTERN "*.html" 
        PATTERN "*.png" 
        PATTERN "*.gif" 
        PATTERN "*.xsl" 
        PATTERN "*.css"
        PATTERN "CMakeFiles" EXCLUDE
        PATTERN "src" EXCLUDE
        )
endif (ALPS_BUILD_SOURCE)


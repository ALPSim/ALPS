set(_alps_tool_source_mean "alea/mean.cpp")
set(_alps_tool_source_variance "alea/variance.cpp")
set(_alps_tool_source_archive "archive/archive.cpp")
set(_alps_tool_source_archive_index "archive/archive_index.cpp")
set(_alps_tool_source_archive_node "archive/archive_node.cpp")
set(_alps_tool_source_archive_plot "archive/archive_plot.cpp")
set(_alps_tool_source_archive_sqlite "archive/archive_sqlite.cpp")
set(_alps_tool_source_archive_xml "archive/archive_xml.cpp")
set(_alps_tool_source_compactrun "scheduler/compactrun.C")
set(_alps_tool_source_convert2xml "scheduler/convert2xml.C")
set(_alps_tool_source_lattice2xml "lattice/lattice2xml.C")
set(_alps_tool_source_maxent "../apps/maxent/cli/maxent.cpp")
set(_alps_tool_source_maxent_helper "../apps/maxent/src/maxent_helper.cpp")
set(_alps_tool_source_maxent_linear_grid_numeric "legacy/maxent_linear_grid_numeric.cpp")
set(_alps_tool_source_maxent_parms "../apps/maxent/src/maxent_parms.cpp")
set(_alps_tool_source_maxent_simulation "../apps/maxent/src/maxent_simulation.cpp")
set(_alps_tool_source_p2h5 "parameters/p2h5.cpp")
set(_alps_tool_source_parameter2hdf5 "parameters/parameter2hdf5.C")
set(_alps_tool_source_parameter2xml "parameters/parameter2xml.C")
set(_alps_tool_source_pconfig "diagnostics/pconfig.C")
set(_alps_tool_source_pevaluate "parapack/pevaluate.C")
set(_alps_tool_source_poutput "parapack/poutput.C")
set(_alps_tool_source_printgraph "lattice/printgraph.C")
set(_alps_tool_source_snap2vtk "scheduler/snap2vtk.C")
set(_alps_tool_source_txt2archive "archive/txt2archive.C")
set(_alps_tool_source_xml2archive "archive/xml2archive.C")
#  Copyright Matthias Troyer and Synge Todo 2009-2013.
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


#
# basic C++ tools
# 


add_subdirectory(alea)

  #
  # Script to setup ALPS_HOME, PATH, LD_LIBRARY_PATH, PYTHONPATH
  #

  if(UNIX AND NOT WIN32)
    configure_file(legacy/alpsvars.sh.in ${CMAKE_CURRENT_BINARY_DIR}/alpsvars.sh @ONLY)
    configure_file(legacy/alpsvars.csh.in ${CMAKE_CURRENT_BINARY_DIR}/alpsvars.csh @ONLY)
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/alpsvars.sh ${CMAKE_CURRENT_BINARY_DIR}/alpsvars.csh DESTINATION bin COMPONENT tools)
  endif(UNIX AND NOT WIN32)


  set (TOOLS parameter2xml parameter2hdf5 convert2xml compactrun lattice2xml pevaluate pconfig
             poutput printgraph txt2archive xml2archive snap2vtk)
  foreach(program ${TOOLS})
    add_executable(${program} ${_alps_tool_source_${program}})
    target_link_libraries(${program} alps)
    install(TARGETS ${program} RUNTIME DESTINATION bin COMPONENT tools) 
  endforeach(program)

  #
  # ngs tools
  # 
  foreach(program p2h5)
    add_executable(${program} ${_alps_tool_source_${program}})
    target_link_libraries(${program} alps)
    install(TARGETS ${program} RUNTIME DESTINATION bin COMPONENT tools) 
  endforeach(program)


if(SQLite_FOUND)
  include_directories(${SQLite_INCLUDE_DIR})
  add_executable(archive archive/archive.cpp archive/archive_index.cpp archive/archive_node.cpp archive/archive_plot.cpp
    archive/archive_sqlite.cpp archive/archive_xml.cpp)
  target_link_libraries(archive alps ${SQLite_LIBRARIES})
  install(TARGETS archive RUNTIME DESTINATION bin  COMPONENT tools)
endif(SQLite_FOUND)

#
# legacy shell scripts that stay around for backwards compatibility until the new solutions are there
#

  if(UNIX AND NOT WIN32)
    install(PROGRAMS legacy/archivecat legacy/convert2html legacy/convert2text legacy/extractgp legacy/extracthtml
            legacy/extracttext legacy/extractxmgr legacy/plot2gp legacy/plot2html
            legacy/plot2text legacy/plot2xmgr legacy/use_local_stylesheet legacy/xslttransform 
            DESTINATION bin COMPONENT tools)
  else(UNIX AND NOT WIN32)
    install(PROGRAMS legacy/msxsl.exe DESTINATION bin COMPONENT tools)
  endif(UNIX AND NOT WIN32)

  #
  # Analytical continuation with MaxEnt
  #

  if(LAPACK_FOUND)
    add_definitions(${LAPACK_DEFINITIONS} -DBIND_FORTRAN_LOWERCASE_UNDERSCORE)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
    add_executable(maxent ../apps/maxent/cli/maxent.cpp ../apps/maxent/src/maxent_helper.cpp ../apps/maxent/src/maxent_simulation.cpp ../apps/maxent/src/maxent_parms.cpp)
    target_link_libraries(maxent alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
    install(TARGETS maxent RUNTIME DESTINATION bin COMPONENT tools)

    # Regression test: linear-FREQUENCY_GRID runs must produce a finite
    # spectrum (guards the t_array_ knot fill in ../apps/maxent/src/maxent_parms.cpp).
    add_executable(maxent_linear_grid_numeric legacy/maxent_linear_grid_numeric.cpp
                   ../apps/maxent/src/maxent_helper.cpp ../apps/maxent/src/maxent_simulation.cpp ../apps/maxent/src/maxent_parms.cpp)
    target_link_libraries(maxent_linear_grid_numeric alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY})
    add_alps_test(maxent_linear_grid_numeric)
   endif(LAPACK_FOUND)

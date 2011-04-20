#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if (ALPS_BUILD_NGS AND ALPS_BUILD_EXAMPLES)

    add_executable(ngs_single main_single.cpp)
    add_dependencies(ngs_single alps)
    target_link_libraries(ngs_single alps)

    add_test(ngs_single
        ${CMAKE_COMMAND}
            -Dsourcedir=${PROJECT_SOURCE_DIR}
            -Dbinarydir=${PROJECT_BINARY_DIR}
            -Dcurrentdir=${CMAKE_CURRENT_SOURCE_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/run_ngs_single_test.cmake
    )

    if (MPI_FOUND)
    
        add_executable(ngs_mpi main_mpi.cpp)
        add_dependencies(ngs_mpi alps)
        target_link_libraries(ngs_mpi alps)
    
        add_test(ngs_mpi
            ${CMAKE_COMMAND}
                -Dsourcedir=${PROJECT_SOURCE_DIR}
                -Dbinarydir=${PROJECT_BINARY_DIR}
                -Dcmddir=${CMAKE_CURRENT_SOURCE_DIR}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/run_ngs_mpi_test.cmake
        )

    endif(MPI_FOUND)
endif (ALPS_BUILD_NGS AND ALPS_BUILD_EXAMPLES)

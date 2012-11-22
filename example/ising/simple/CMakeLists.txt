#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if (ALPS_BUILD_EXAMPLES)

    add_executable(example_ising_single single)
    add_dependencies(example_ising_single alps)
    target_link_libraries(example_ising_single alps)

    add_test(example_ising_single
        ${CMAKE_COMMAND}
            -Dcmd=ngs_ising_single
            -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
            -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
            -Dbindir=${PROJECT_BINARY_DIR}/bin
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ngs_ising_test.cmake
    )

    if (MPI_FOUND)

        add_executable(example_ising_mpi mpi.cpp)
        add_dependencies(example_ising_mpi alps)
        target_link_libraries(example_ising_mpi alps)

        add_test(ngs_ising_mpi
            ${CMAKE_COMMAND}
                -Dcmd=ngs_ising_mpi
                -Dmpiexec=${MPIEXEC}
                -Dmpiexec_numproc_flag=${MPIEXEC_NUMPROC_FLAG}
                -Dprocs=4
                -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
                -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
                -Dbindir=${PROJECT_BINARY_DIR}/bin
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ngs_ising_test.cmake
        )

    endif(MPI_FOUND)

endif (ALPS_BUILD_EXAMPLES)

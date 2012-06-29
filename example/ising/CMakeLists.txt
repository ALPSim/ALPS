#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if (ALPS_BUILD_NGS AND ALPS_BUILD_EXAMPLES)

    add_executable(ngs_ising_single main_single.cpp)
    add_dependencies(ngs_ising_single alps)
    target_link_libraries(ngs_ising_single alps)

    add_test(ngs_ising_single
        ${CMAKE_COMMAND}
            -Dsourcedir=${PROJECT_SOURCE_DIR}
            -Dbinarydir=${PROJECT_BINARY_DIR}
            -Dcurrentdir=${CMAKE_CURRENT_SOURCE_DIR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ngs_ising_single_test.cmake
    )

# TODO: implement
#    add_executable(ngs_ising_threaded main_threaded.cpp)
#    add_dependencies(ngs_ising_threaded alps)
#    target_link_libraries(ngs_ising_threaded alps)

#    add_test(ngs_ising_threaded
#        ${CMAKE_COMMAND}
#            -Dsourcedir=${PROJECT_SOURCE_DIR}
#            -Dbinarydir=${PROJECT_BINARY_DIR}
#            -Dcurrentdir=${CMAKE_CURRENT_SOURCE_DIR}
#            -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ngs_ising_threaded_test.cmake
#    )

    if (MPI_FOUND)
    
        add_executable(ngs_ising_mpi main_mpi.cpp)
        add_dependencies(ngs_ising_mpi alps)
        target_link_libraries(ngs_ising_mpi alps)
    
        add_test(ngs_ising_mpi
            ${CMAKE_COMMAND}
                -Dsourcedir=${PROJECT_SOURCE_DIR}
                -Dbinarydir=${PROJECT_BINARY_DIR}
                -Dcurrentdir=${CMAKE_CURRENT_SOURCE_DIR}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ngs_ising_mpi_test.cmake
        )

# TODO: Implement
#        add_executable(ngs_ising_mpi_threaded main_mpi_threaded.cpp)
#        add_dependencies(ngs_ising_mpi_threaded alps)
#        target_link_libraries(ngs_ising_mpi_threaded alps)
    
#        add_test(ngs_ising_mpi_threaded
#            ${CMAKE_COMMAND}
#                -Dsourcedir=${PROJECT_SOURCE_DIR}
#                -Dbinarydir=${PROJECT_BINARY_DIR}
#                -Dcurrentdir=${CMAKE_CURRENT_SOURCE_DIR}
#                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ngs_ising_mpi_threaded_test.cmake
#        )

        add_executable(ngs_ising_mpi_pscan main_mpi_pscan.cpp)
        add_dependencies(ngs_ising_mpi_pscan alps)
        target_link_libraries(ngs_ising_mpi_pscan alps)
    
        add_test(ngs_ising_mpi_pscan
            ${CMAKE_COMMAND}
                -Dsourcedir=${PROJECT_SOURCE_DIR}
                -Dbinarydir=${PROJECT_BINARY_DIR}
                -Dcurrentdir=${CMAKE_CURRENT_SOURCE_DIR}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ngs_ising_mpi_pscan_test.cmake
        )

        add_executable(ngs_ising_mpi2 main_mpi2.cpp)
        add_dependencies(ngs_ising_mpi2 alps)
        target_link_libraries(ngs_ising_mpi2 alps)
   
        add_test(ngs_ising_mpi2
            ${CMAKE_COMMAND}
                -Dsourcedir=${PROJECT_SOURCE_DIR}
                -Dbinarydir=${PROJECT_BINARY_DIR}
                -Dcurrentdir=${CMAKE_CURRENT_SOURCE_DIR}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ngs_ising_mpi_test.cmake
        )

    endif(MPI_FOUND)

    if (ALPS_HAVE_PYTHON)

        set(ALPS_SHARED_CPPFLAGS ALPS_EXPORTS=1)
        set(ALPS_STATIC_CPPFLAGS "")

        python_add_module(pyising pyising.cpp)

        if(BUILD_SHARED_LIBS)
            set_target_properties(pyising PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
            set_target_properties(pyising PROPERTIES PREFIX "")
            if(WIN32 AND NOT UNIX)
                set_target_properties(pyising PROPERTIES SUFFIX ".pyd")
            endif(WIN32 AND NOT UNIX)
        endif(BUILD_SHARED_LIBS)
        target_link_libraries(pyising alps ${LINK_LIBRARIES})

#        install(TARGETS pyising
#                COMPONENT python
#                RUNTIME DESTINATION example/ising
#                ARCHIVE DESTINATION example/ising
#                LIBRARY DESTINATION example/ising)

        add_test(ngs_ising_python_single
            ${CMAKE_COMMAND}
                -Dpython_interpreter=${PYTHON_INTERPRETER} 
                -Dsourcedir=${PROJECT_SOURCE_DIR}
                -Dbinarydir=${PROJECT_BINARY_DIR}
                -Dcurrentdir=${CMAKE_CURRENT_SOURCE_DIR}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ngs_ising_python_single_test.cmake
        )

        add_test(ngs_ising_python_native
            ${CMAKE_COMMAND}
                -Dpython_interpreter=${PYTHON_INTERPRETER}
                -Dsourcedir=${PROJECT_SOURCE_DIR}
                -Dbinarydir=${PROJECT_BINARY_DIR}
                -Dcurrentdir=${CMAKE_CURRENT_SOURCE_DIR}
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ngs_ising_python_native_test.cmake
        )

        if(MPI_FOUND)
            add_test(ngs_ising_python_native_mpi
                ${CMAKE_COMMAND}
                    -Dpython_interpreter=${PYTHON_INTERPRETER}
                    -Dsourcedir=${PROJECT_SOURCE_DIR}
                    -Dbinarydir=${PROJECT_BINARY_DIR}
                    -Dcurrentdir=${CMAKE_CURRENT_SOURCE_DIR}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ngs_ising_python_native_mpi_test.cmake
            )
        endif(MPI_FOUND)

    endif (ALPS_HAVE_PYTHON)

endif (ALPS_BUILD_NGS AND ALPS_BUILD_EXAMPLES)

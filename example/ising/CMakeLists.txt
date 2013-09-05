#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

# TODO: fix it ...
if (FALSE AND ALPS_BUILD_EXAMPLES)

    add_executable(alpscontrol alpscontrol.cpp)
    add_dependencies(alpscontrol alps)
    target_link_libraries(alpscontrol alps)

    foreach (name single controlthread_tiny server_tiny)
    
        add_executable(ngs_ising_${name} ${name}.cpp)
        add_dependencies(ngs_ising_${name} alps)
        target_link_libraries(ngs_ising_${name} alps)

        add_test(ngs_ising_${name}
            ${CMAKE_COMMAND}
                -Dcmd=ngs_ising_${name}
                -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
                -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
                -Dbindir=${PROJECT_BINARY_DIR}/bin
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ngs_ising_test.cmake
        )

    endforeach(name)

    foreach (name singlexml)
    
        add_executable(ngs_ising_${name} ${name}.cpp)
        add_dependencies(ngs_ising_${name} alps)
        target_link_libraries(ngs_ising_${name} alps)

        add_test(ngs_ising_${name}
            ${CMAKE_COMMAND}
                -Dcmd=ngs_ising_${name}
                -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
                -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
                -Dbindir=${PROJECT_BINARY_DIR}/bin
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ngs_ising_test_xml.cmake
        )

    endforeach(name)

    add_executable(ngs_parapack_ising parapack.cpp)
    # add_dependencies(ngs_parapack_ising alps)
    target_link_libraries(ngs_parapack_ising alps)

    if (MPI_FOUND)

        SET(ALPS_NGS_EXAMPLES_MPI_EXEC mpi mpi_pscan mpi_controlthread_tiny mpi_server_tiny parallel)
        if (ALPS_NGS_OPENMPI_ULFM)
            SET(ALPS_NGS_EXAMPLES_MPI_EXEC ${ALPS_NGS_EXAMPLES_MPI_EXEC} mpi_ulfm)
        endif(ALPS_NGS_OPENMPI_ULFM)

        foreach (name ${ALPS_NGS_EXAMPLES_MPI_EXEC})

        
            add_executable(ngs_ising_${name} ${name}.cpp)
            add_dependencies(ngs_ising_${name} alps)
            target_link_libraries(ngs_ising_${name} alps)

            add_test(ngs_ising_${name}
                ${CMAKE_COMMAND}
                    -Dcmd=ngs_ising_${name}
                    -Dmpiexec=${MPIEXEC}
                    -Dmpiexec_numproc_flag=${MPIEXEC_NUMPROC_FLAG}
                    -Dprocs=4
                    -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
                    -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
                    -Dbindir=${PROJECT_BINARY_DIR}/bin
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ngs_ising_test.cmake
            )

        endforeach(name)

    endif(MPI_FOUND)

# TODO: use mcbase_ng and not mcbase_ng
    if (ALPS_HAVE_PYTHON AND DISABLE_SINCE_MCBASE_NG_IS_NOT_EXPORTED)

        set(ALPS_SHARED_CPPFLAGS ALPS_EXPORTS=1)
        set(ALPS_STATIC_CPPFLAGS "")

        python_add_module(ngsising_c src/pyising.cpp)

        if(BUILD_SHARED_LIBS)
            set_target_properties(ngsising_c PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
            set_target_properties(ngsising_c PROPERTIES PREFIX "")
            if(WIN32 AND NOT UNIX)
                set_target_properties(ngsising_c PROPERTIES SUFFIX ".pyd")
            endif(WIN32 AND NOT UNIX)
        endif(BUILD_SHARED_LIBS)
        target_link_libraries(ngsising_c alps ${LINK_LIBRARIES})

        install(TARGETS ngsising_c
                COMPONENT python
                RUNTIME DESTINATION example/ising
                ARCHIVE DESTINATION example/ising
                LIBRARY DESTINATION example/ising)

        foreach (name export native)

            add_test(ngs_ising_python_${name}
                ${CMAKE_COMMAND}
                    -Dcmd=${name}
                    -Dpython_interpreter=${PYTHON_INTERPRETER}
                    -Dsourcedir=${PROJECT_SOURCE_DIR}
                    -Dbinarydir=${PROJECT_BINARY_DIR}
                    -Dbindir=${PROJECT_BINARY_DIR}/bin
                    -Dcurrentdir=${CMAKE_CURRENT_SOURCE_DIR}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ngs_ising_python_test.cmake
            )
            
        endforeach(name)

        if(MPI_FOUND)

            add_test(ngs_ising_python_mpi_native
                ${CMAKE_COMMAND}
                    -Dcmd=native
                    -Dpython_interpreter=${PYTHON_INTERPRETER}
                    -Dmpiexec=${MPIEXEC}
                    -Dmpiexec_numproc_flag=${MPIEXEC_NUMPROC_FLAG}
                    -Dprocs=4
                    -Dsourcedir=${PROJECT_SOURCE_DIR}
                    -Dbinarydir=${PROJECT_BINARY_DIR}
                    -Dbindir=${PROJECT_BINARY_DIR}/bin
                    -Dcurrentdir=${CMAKE_CURRENT_SOURCE_DIR}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_ngs_ising_python_test.cmake
            )

        endif(MPI_FOUND)

    endif (ALPS_HAVE_PYTHON AND DISABLE_SINCE_MCBASE_NG_IS_NOT_EXPORTED)

    add_subdirectory(simple)

endif (FALSE AND ALPS_BUILD_EXAMPLES)

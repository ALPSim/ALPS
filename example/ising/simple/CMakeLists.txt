#  Copyright Matthias Troyer 2009 - 2010.
#  Distributed under the Boost Software License, Version 1.0.
#      (See accompanying file LICENSE_1_0.txt or copy at
#          http://www.boost.org/LICENSE_1_0.txt)

if (ALPS_BUILD_EXAMPLES)


    set(ALPS_EXAMPLES_NGS_SIMPLE single intobs)

    if(ALPS_NGS_USE_NEW_ALEA)
    
    set(ALPS_EXAMPLES_NGS_SIMPLE ${ALPS_EXAMPLES_NGS_SIMPLE} multi_array)

    endif(ALPS_NGS_USE_NEW_ALEA)
    
    foreach (name ${ALPS_EXAMPLES_NGS_SIMPLE})
    
        add_executable(example_ising_${name} ${name}.cpp)
        add_dependencies(example_ising_${name} alps)
        target_link_libraries(example_ising_${name} alps)

        add_test(example_ising_${name}
            ${CMAKE_COMMAND}
                -Dcmd=example_ising_${name}
                -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
                -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
                -Dbindir=${PROJECT_BINARY_DIR}/bin
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_test.cmake
        )

    endforeach(name)

    if (MPI_FOUND)

        add_executable(example_ising_mpi mpi.cpp)
        add_dependencies(example_ising_mpi alps)
        target_link_libraries(example_ising_mpi alps)

        add_test(example_ising_mpi
            ${CMAKE_COMMAND}
                -Dcmd=example_ising_mpi
                -Dmpiexec=${MPIEXEC}
                -Dmpiexec_numproc_flag=${MPIEXEC_NUMPROC_FLAG}
                -Dprocs=4
                -Dsourcedir=${CMAKE_CURRENT_SOURCE_DIR}
                -Dbinarydir=${CMAKE_CURRENT_BINARY_DIR}
                -Dbindir=${PROJECT_BINARY_DIR}/bin
                -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_test.cmake
        )

    endif(MPI_FOUND)

    if (ALPS_HAVE_PYTHON)

        set(ALPS_SHARED_CPPFLAGS ALPS_EXPORTS=1)
        set(ALPS_STATIC_CPPFLAGS "")

        python_add_module(exampleising_c export2py.cpp)

        if(BUILD_SHARED_LIBS)
            set_target_properties(exampleising_c PROPERTIES COMPILE_DEFINITIONS "${ALPS_SHARED_CPPFLAGS}")
            set_target_properties(exampleising_c PROPERTIES PREFIX "")
            if(WIN32 AND NOT UNIX)
                set_target_properties(exampleising_c PROPERTIES SUFFIX ".pyd")
            endif(WIN32 AND NOT UNIX)
        endif(BUILD_SHARED_LIBS)
        target_link_libraries(exampleising_c alps ${LINK_LIBRARIES})

        install(TARGETS exampleising_c
                COMPONENT python
                RUNTIME DESTINATION example/ising/simple
                ARCHIVE DESTINATION example/ising/simple
                LIBRARY DESTINATION example/ising/simple)

        foreach (name export native)

            add_test(example_ising_python_${name}
                ${CMAKE_COMMAND}
                    -Dcmd=${name}
                    -Dpython_interpreter=${PYTHON_INTERPRETER}
                    -Dsourcedir=${PROJECT_SOURCE_DIR}
                    -Dbinarydir=${PROJECT_BINARY_DIR}
                    -Dbindir=${PROJECT_BINARY_DIR}/bin
                    -Dcurrentdir=${CMAKE_CURRENT_SOURCE_DIR}
                    -P ${CMAKE_CURRENT_SOURCE_DIR}/tests/run_python_test.cmake
            )

        endforeach (name)

    endif(ALPS_HAVE_PYTHON)

endif (ALPS_BUILD_EXAMPLES)

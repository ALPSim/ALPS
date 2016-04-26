if(NOT LAPACK_FOUND)
  message(STATUS "MPS will not be built since lapack library is not found")
else(NOT LAPACK_FOUND)
  ### DEFINITIONS AND INCLUDES VARIABLES (ready to be exported)
  set(ALPS_MPS_DEFINITIONS)
  set(ALPS_MPS_INCLUDE_DIRS)
  set(ALPS_MPS_LINKER_FLAGS)
  
  ### SYMMETRY MANAGEMENT
  include(symmetries.cmake)
  set_symmetry_difinitions(SYMMETRY_DEFINITIONS)
  list(APPEND ALPS_MPS_DEFINITIONS ${SYMMETRY_DEFINITIONS})

  ### PREPARE SETTINGS
  list(APPEND ALPS_MPS_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR})
  list(APPEND ALPS_MPS_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/framework)
  list(APPEND ALPS_MPS_INCLUDE_DIRS ${CMAKE_CURRENT_BINARY_DIR}/framework)
  
  list(APPEND ALPS_MPS_DEFINITIONS "-DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING")
  list(APPEND ALPS_MPS_DEFINITIONS "-DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING")
  
  set(ALPS_MPS_LINKER_FLAGS "${LAPACK_LINKER_FLAGS}")  

  if(OPENMP_FOUND)
    list(APPEND ALPS_MPS_DEFINITIONS "-DMAQUIS_OPENMP")
  endif(OPENMP_FOUND)
  

  ## CLEANUP
  ## clean binaries created in previous versions so that config dirs can be created
  if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/mps_optim" AND NOT IS_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mps_optim")
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/mps_optim")
  endif()
  if(EXISTS "${CMAKE_CURRENT_BINARY_DIR}/mps_evolve" AND NOT IS_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/mps_evolve")
    file(REMOVE "${CMAKE_CURRENT_BINARY_DIR}/mps_evolve")
  endif()


  ### FRAMEWORK
  set(FRAMEWORK_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/framework/dmrg)
  set(FRAMEWORK_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/framework/dmrg)

  configure_file( ${FRAMEWORK_SOURCE_DIR}/block_matrix/symmetry/nu1.h.in
                  ${FRAMEWORK_BINARY_DIR}/block_matrix/symmetry/nu1.h )

  set(MPS_UTILS_SOURCES ${FRAMEWORK_SOURCE_DIR}/block_matrix/symmetry.cpp
                        ${FRAMEWORK_SOURCE_DIR}/utils/utils.cpp
                        ${FRAMEWORK_SOURCE_DIR}/utils/time_stopper.cpp)

  set(MPS_MODELS_SOURCES "${FRAMEWORK_SOURCE_DIR}/models/factory/lattice_factory.cpp")
  get_symm_files(APPEND MPS_MODELS_SOURCES "${FRAMEWORK_SOURCE_DIR}/models/factory/model_factory_{SYMM}.cpp")
  get_symm_files(APPEND MPS_MODELS_SOURCES "${FRAMEWORK_SOURCE_DIR}/models/alps/symm_handler_impl/model_{SYMM}.cpp")
  
  ### SETTINGS
  add_definitions(${ALPS_MPS_DEFINITIONS})
  include_directories(${ALPS_MPS_INCLUDE_DIRS})
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${ALPS_MPS_LINKER_FLAGS}")
  set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ALPS_MPS_LINKER_FLAGS}")

  
  
  ### LIBRARY TARGETS
  add_library(mps_utils SHARED ${MPS_UTILS_SOURCES})
  target_link_libraries(mps_utils alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${HDF5_LIBRARIES})
  add_library(mps_models SHARED ${MPS_MODELS_SOURCES})
  target_link_libraries(mps_models mps_utils alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${HDF5_LIBRARIES})
  #### SCHEDULER
  add_library(mps_pscan STATIC libpscan/scheduler.cpp libpscan/options.cpp)
  
  ### INSTALL LIBS
  install(TARGETS mps_models mps_utils DESTINATION lib EXPORT ALPSMPS-targets)
  
  ### APPLICATIONS
  # ${ALPS_BOOST_LIBRARY_NAME} ${HDF5_LIBRARIES} included in `alps`
  set(MYAPP_LIBRARIES mps_pscan mps_models mps_utils alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${HDF5_LIBRARIES})
  
  add_subdirectory(mps_eigenstate_sims)
  add_subdirectory(mps_evolve)
  add_subdirectory(mps_overlap)
  add_subdirectory(examples)
  
  
  if(ALPS_INSTALL_MPS_FRAMEWORK)
    ## Install headers and libs
    install(DIRECTORY framework/dmrg framework/utils
            DESTINATION include
            COMPONENT headers
            PATTERN ".svn" EXCLUDE
            PATTERN ".DS_Store" EXCLUDE
            PATTERN "*.cpp" EXCLUDE
            PATTERN "*.C" EXCLUDE
            PATTERN "*.c" EXCLUDE
            PATTERN "*.cu" EXCLUDE
            PATTERN "*.in" EXCLUDE
            PATTERN "CMakeLists.txt" EXCLUDE
            PATTERN "Makefile*" EXCLUDE
            )
    install(FILES ${FRAMEWORK_BINARY_DIR}/block_matrix/symmetry/nu1.h
            DESTINATION include/dmrg/block_matrix/symmetry COMPONENT headers)
    
    ## Install ALPSMPSConfig.cmake and targets
    set(CONF_INCLUDE_DIRS ${CMAKE_INSTALL_PREFIX}/include) # TODO: add hdf5&others explicitely
    set(CONF_LIBRARY_DIRS ${CMAKE_INSTALL_PREFIX}/lib) # TODO: add hdf5&others explicitely
    set(CONF_SYMMETRY_CMAKE ${CMAKE_INSTALL_PREFIX}/share/alps/symmetries.cmake)
    set(CONF_LIBRARIES mps_models mps_utils alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${HDF5_LIBRARIES})
    
    configure_file(ALPS_MPSConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/ALPS_MPSConfig.cmake @ONLY)
    
    install(FILES symmetries.cmake DESTINATION share/alps)
    install(FILES ${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/ALPS_MPSConfig.cmake DESTINATION share/alps)
    # TODO: Export ALPS targets that are required by this.
    # install(EXPORT ALPSMPS-targets DESTINATION share/alps)
  endif(ALPS_INSTALL_MPS_FRAMEWORK)
  
endif(NOT LAPACK_FOUND)

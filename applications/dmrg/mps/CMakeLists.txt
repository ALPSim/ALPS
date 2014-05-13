if(NOT LAPACK_FOUND)
  message(STATUS "MPS will not be built since lapack library is not found")
else(NOT LAPACK_FOUND)
  ### SYMMETRY MANAGEMENT
  include(symmetries.cmake)
  set_symmetry_difinitions(SYMMETRY_DEFINITIONS)
  add_definitions(${SYMMETRY_DEFINITIONS})

  ### SETTINGS
  include_directories(${CMAKE_CURRENT_SOURCE_DIR})
  include_directories(${CMAKE_CURRENT_SOURCE_DIR}/framework)
  include_directories(${CMAKE_CURRENT_BINARY_DIR}/framework)

  add_definitions(-DDISABLE_MATRIX_ELEMENT_ITERATOR_WARNING -DALPS_DISABLE_MATRIX_ELEMENT_ITERATOR_WARNING)

  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")

  if(OPENMP_FOUND)
    add_definitions(-DMAQUIS_OPENMP)
  endif(OPENMP_FOUND)

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

  add_library(mps_utils STATIC ${MPS_UTILS_SOURCES})
  add_library(mps_models STATIC ${MPS_MODELS_SOURCES})


  #### SCHEDULER
  add_library(mps_pscan STATIC libpscan/scheduler.cpp libpscan/options.cpp)


  ### APPLICATIONS
  # ${ALPS_BOOST_LIBRARY_NAME} ${HDF5_LIBRARIES} included in `alps`
  set(MYAPP_LIBRARIES mps_pscan mps_models mps_utils alps ${LAPACK_LIBRARY} ${BLAS_LIBRARY} ${HDF5_LIBRARIES})
  
  add_subdirectory(mps_optim)
  add_subdirectory(mps_evolve)

endif(NOT LAPACK_FOUND)

# Scheduler source ownership.
target_sources(alps PRIVATE
  src/scheduler/abstract_task.C
  src/scheduler/convertxdr.C
  src/scheduler/factory.C
  src/scheduler/info.C
  src/scheduler/master_scheduler.C
  src/scheduler/measurement_operators.C
  src/scheduler/montecarlo.C
  src/scheduler/mpp_scheduler.C
  src/scheduler/options.C
  src/scheduler/remote_task.C
  src/scheduler/remote_worker.C
  src/scheduler/scheduler.C
  src/scheduler/serial_scheduler.C
  src/scheduler/signal.C
  src/scheduler/single_scheduler.C
  src/scheduler/slave_task.C
  src/scheduler/task.C
  src/scheduler/worker.C
  src/scheduler/workertask.C)

if(ALPS_NGS_OPENMPI_ULFM)
  if(NOT ALPS_ENABLE_MPI)
    message(FATAL_ERROR "ALPS_NGS_OPENMPI_ULFM requires ALPS_ENABLE_MPI=ON")
  endif()
  target_sources(alps PRIVATE src/ngs/ulfm.cpp)
endif()

target_sources(alps_headers INTERFACE
  FILE_SET scheduler_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/ngs/detail/tcpsession.hpp
    include/alps/ngs/scheduler/proto/controlthreadsim.hpp
    include/alps/ngs/scheduler/proto/mcbase.hpp
    include/alps/ngs/scheduler/proto/mpisim.hpp
    include/alps/ngs/scheduler/proto/mpisim_ulfm.hpp
    include/alps/ngs/scheduler/proto/tcpserver.hpp
    include/alps/ngs/scheduler/tcpserver.hpp
    include/alps/ngs/ulfm.hpp
    include/alps/scheduler.h
    include/alps/scheduler/convert.h
    include/alps/scheduler/diag.hpp
    include/alps/scheduler/factory.h
    include/alps/scheduler/info.h
    include/alps/scheduler/measurement_operators.h
    include/alps/scheduler/montecarlo.h
    include/alps/scheduler/options.h
    include/alps/scheduler/scheduler.h
    include/alps/scheduler/signal.hpp
    include/alps/scheduler/task.h
    include/alps/scheduler/types.h
    include/alps/scheduler/worker.h)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

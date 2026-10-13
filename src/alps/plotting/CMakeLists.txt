# Plot output combines XML with the older Parameters data model.
target_sources(alps_headers INTERFACE
  FILE_SET plotting_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/plot.h)

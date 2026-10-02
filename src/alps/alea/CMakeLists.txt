# Alea source ownership.
target_sources(alps PRIVATE
  src/alea/nan.C
  src/alea/observable.C
  src/alea/observablefactory.C
  src/alea/observableset.C
  src/ngs/mcobservable.cpp
  src/ngs/mcobservables.cpp
  src/ngs/mcresult.cpp
  src/ngs/mcresults.cpp
  src/ngs/observablewrappers.cpp)

target_sources(alps_headers INTERFACE
  FILE_SET alea_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/alea.h
    include/alps/alea/abstractbinning.h
    include/alps/alea/abstractsimpleobservable.h
    include/alps/alea/abstractsimpleobservable.ipp
    include/alps/alea/convergence.hpp
    include/alps/alea/detailedbinning.h
    include/alps/alea/histogram.h
    include/alps/alea/histogramdata.h
    include/alps/alea/histogrameval.h
    include/alps/alea/mcanalyze.hpp
    include/alps/alea/mcdata.hpp
    include/alps/alea/nan.h
    include/alps/alea/nobinning.h
    include/alps/alea/observable.h
    include/alps/alea/observable_fwd.hpp
    include/alps/alea/observablefactory.h
    include/alps/alea/observableset.h
    include/alps/alea/observableset_p.h
    include/alps/alea/output_helper.h
    include/alps/alea/recordableobservable.h
    include/alps/alea/signedobservable.h
    include/alps/alea/simplebinning.h
    include/alps/alea/simpleobsdata.h
    include/alps/alea/simpleobservable.h
    include/alps/alea/simpleobservable.ipp
    include/alps/alea/simpleobseval.h
    include/alps/alea/simpleobseval.ipp
    include/alps/alea/type_tag.hpp
    include/alps/alea/value_with_error.hpp
    include/alps/ngs/lib/mcresult_impl_base.ipp
    include/alps/ngs/lib/mcresult_impl_derived.ipp
    include/alps/ngs/mcobservable.hpp
    include/alps/ngs/mcobservables.hpp
    include/alps/ngs/mcresult.hpp
    include/alps/ngs/mcresults.hpp
    include/alps/ngs/observablewrappers.hpp)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

alps_register_module(NAME alea HEADER_SET alea_headers
  TARGETS alps
  DEPENDS configuration containers hdf5 numerics osiris utilities xml
)

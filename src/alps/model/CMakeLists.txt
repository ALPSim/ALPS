# Model source ownership.
target_sources(alps PRIVATE
  src/model/bondterm.C
  src/model/globaloperator.C
  src/model/modellibrary.C
  src/model/siteterm.C)

target_sources(alps_headers INTERFACE
  FILE_SET model_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/model.h
    include/alps/model/basisdescriptor.h
    include/alps/model/basisstates.h
    include/alps/model/blochbasisstates.h
    include/alps/model/bondoperator.h
    include/alps/model/bondterm.h
    include/alps/model/default_term.h
    include/alps/model/globaloperator.h
    include/alps/model/half_integer.h
    include/alps/model/hamiltonian.h
    include/alps/model/hamiltonian_matrix.hpp
    include/alps/model/integer_state.h
    include/alps/model/model_helper.h
    include/alps/model/modellibrary.h
    include/alps/model/operator.h
    include/alps/model/operatordescriptor.h
    include/alps/model/operatorsubstitution.h
    include/alps/model/quantumnumber.h
    include/alps/model/sign.h
    include/alps/model/sitebasisdescriptor.h
    include/alps/model/sitebasisstates.h
    include/alps/model/siteoperator.h
    include/alps/model/sitestate.h
    include/alps/model/siteterm.h
    include/alps/model/substitute.h)

if(ALPS_BUILD_TESTING)
  add_subdirectory(tests)
endif()

alps_register_module(NAME model HEADER_SET model_headers
  TARGETS alps
  DEPENDS configuration containers expression lattice legacy_parameters numerics utilities xml
)

# Expression source ownership.
target_sources(alps PRIVATE
  src/expression/evaluator.C)

target_sources(alps_headers INTERFACE
  FILE_SET expression_headers TYPE HEADERS
  BASE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/include"
  FILES
    include/alps/expression.h
    include/alps/expression/block.h
    include/alps/expression/evaluatable.h
    include/alps/expression/evaluate.h
    include/alps/expression/evaluate_helper.h
    include/alps/expression/evaluator.h
    include/alps/expression/expression.h
    include/alps/expression/expression_fwd.h
    include/alps/expression/factor.h
    include/alps/expression/function.h
    include/alps/expression/number.h
    include/alps/expression/parameterevaluator.h
    include/alps/expression/symbol.h
    include/alps/expression/term.h
    include/alps/expression/traits.h)

alps_register_module(NAME expression HEADER_SET expression_headers
  TARGETS alps
  DEPENDS configuration legacy_parameters numerics random utilities xml
)

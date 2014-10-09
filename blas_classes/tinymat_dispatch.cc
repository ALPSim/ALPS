//by Emanuel
//construct functions with name matmult_p_q_r(C, A, B) where C has size p,q; A
//has size p,r and B has size r,q.
//use a good compiler with unrolling and sse2/sse3 to get fast execution for
//small matrices here. If possible use IPO to inline these functions later.
#include <iostream>
//this fixes an Apple bug in tvmet
#undef __APPLE_CC__
//replacement for blitz: tiny matrix and vector library
#include <tvmet/Matrix.h>
#include <tvmet/Vector.h>

using namespace tvmet;
//boost PP headers
#include <boost/preprocessor/seq/elem.hpp>
#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/cat.hpp>

#define MSIZE_P (1)(2)(4)(6)(8)(10)(12)
#define MSIZE_Q (1)(2)(4)(6)(8)(10)(12)
#define MSIZE_R (1)(2)(4)(6)(8)(10)(12)

#define MLOOP(unused, product) \
  BOOST_PP_CAT(\
    BOOST_PP_CAT(\
      BOOST_PP_CAT(\
        BOOST_PP_CAT(\
          BOOST_PP_CAT(\
            void matmult_, BOOST_PP_SEQ_ELEM(0, product)\
          )\
          , _\
        )\
        , BOOST_PP_SEQ_ELEM(1, product)\
      )\
      , _\
    )\
    , BOOST_PP_SEQ_ELEM(2, product)\
  )(const double * A, const double * B, double * C){\
    *((tvmet::Matrix<double, BOOST_PP_SEQ_ELEM(0, product), BOOST_PP_SEQ_ELEM(1, product)> *)C)=\
      prod(*((const tvmet::Matrix<double, BOOST_PP_SEQ_ELEM(0, product), BOOST_PP_SEQ_ELEM(2, product)>*)(A)),\
           *((const tvmet::Matrix<double, BOOST_PP_SEQ_ELEM(2, product), BOOST_PP_SEQ_ELEM(1, product)>*)(B))); \
  }

//that's where we're actually producing the function bodies. To look at it,
//just preprocess and inspect visually ;-)
BOOST_PP_SEQ_FOR_EACH_PRODUCT(MLOOP, (MSIZE_P)(MSIZE_Q)(MSIZE_R))

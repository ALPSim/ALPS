/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2010 - 2011 by Lukas Gamper <gamperl@gmail.com>                   *
 *                                                                                 *
 * Permission is hereby granted, free of charge, to any person obtaining           *
 * a copy of this software and associated documentation files (the “Software”),    *
 * to deal in the Software without restriction, including without limitation       *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,        *
 * and/or sell copies of the Software, and to permit persons to whom the           *
 * Software is furnished to do so, subject to the following conditions:            *
 *                                                                                 *
 * The above copyright notice and this permission notice shall be included         *
 * in all copies or substantial portions of the Software.                          *
 *                                                                                 *
 * THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS         *
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,     *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     *
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          *
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING         *
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER             *
 * DEALINGS IN THE SOFTWARE.                                                       *
 *                                                                                 *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef ALPS_NGS_ALEA_ACCUMULATOR_HEADER
#define ALPS_NGS_ALEA_ACCUMULATOR_HEADER

#include <alps/ngs/alea/result.hpp>
#include <alps/ngs/alea/features.hpp>
#include <alps/ngs/alea/accumulator/arguments.hpp>
#include <alps/ngs/alea/accumulator/accumulator_impl.hpp>

#include <boost/parameter.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits.hpp>

// = = = = N A M E D   P A R A M E T E R   C T O R   D E F I N I T I O N = = = =

namespace alps {
    namespace accumulator  {
        template<
              typename vt  = double //value_type
            , typename features_input = features<tag::mean, tag::error>
            , typename wvt = detail::no_weight_value_type //default: there is no wvt, therefor it's not a weighted accumulator
        >
        class accumulator: public detail::accumulator_impl<
              type_holder<vt, wvt>
            , typename features_input::T0
            , typename features_input::T1
            , typename features_input::T2
            , typename features_input::T3
            , typename features_input::T4
            , typename features_input::T5
            , typename features_input::T6
            , typename features_input::T7
            , typename features_input::T8
            , typename features_input::T9
            , typename features_input::T10
            , typename features_input::T11
            , typename features_input::T12
            , typename features_input::T13
            , typename features_input::T14
            , typename features_input::T15
            , typename features_input::T16
            , typename features_input::T17
            , typename features_input::T18
            , typename features_input::T19
            , typename boost::mpl::if_c<boost::is_same<wvt, detail::no_weight_value_type>::value, void , tag::detail::weight>::type
        > {
            typedef accumulator<vt, features_input, wvt> self_type;

            typedef detail::accumulator_impl<
                  type_holder<vt, wvt>
                , typename features_input::T0
                , typename features_input::T1
                , typename features_input::T2
                , typename features_input::T3
                , typename features_input::T4
                , typename features_input::T5
                , typename features_input::T6
                , typename features_input::T7
                , typename features_input::T8
                , typename features_input::T9
                , typename features_input::T10
                , typename features_input::T11
                , typename features_input::T12
                , typename features_input::T13
                , typename features_input::T14
                , typename features_input::T15
                , typename features_input::T16
                , typename features_input::T17
                , typename features_input::T18
                , typename features_input::T19
                , typename boost::mpl::if_c<boost::is_same<wvt, detail::no_weight_value_type>::value, void , tag::detail::weight>::type
            > base_type;
            
            public:

                typedef result<
                      type_holder<vt, wvt>
                    , typename features_input::T0
                    , typename features_input::T1
                    , typename features_input::T2
                    , typename features_input::T3
                    , typename features_input::T4
                    , typename features_input::T5
                    , typename features_input::T6
                    , typename features_input::T7
                    , typename features_input::T8
                    , typename features_input::T9
                    , typename features_input::T10
                    , typename features_input::T11
                    , typename features_input::T12
                    , typename features_input::T13
                    , typename features_input::T14
                    , typename features_input::T15
                    , typename features_input::T16
                    , typename features_input::T17
                    , typename features_input::T18
                    , typename features_input::T19
                    , typename boost::mpl::if_c<boost::is_same<wvt, detail::no_weight_value_type>::value, void , tag::detail::weight>::type
                > result_type;

                accumulator(accumulator const & arg): base_type(static_cast<base_type const &>(arg)) {}
            
                BOOST_PARAMETER_CONSTRUCTOR(
                accumulator, 
                (base_type),
                keywords,
                    (optional 
                        (_bin_size, *)
                        (_bin_num, *)
                        (_weight_ref, *)
                    )
                )
        };
    }
}
#endif

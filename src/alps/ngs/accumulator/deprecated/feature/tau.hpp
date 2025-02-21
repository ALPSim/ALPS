/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *                                                                                 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations                  *
 *                                                                                 *
 * ALPS Libraries                                                                  *
 *                                                                                 *
 * Copyright (C) 2011 - 2012 by Mario Koenz <mkoenz@ethz.ch>                       *
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

#ifndef ALPS_NGS_ALEA_DETAIL_TAU_IMPLEMENTATION_HEADER
#define ALPS_NGS_ALEA_DETAIL_TAU_IMPLEMENTATION_HEADER

#include <alps/ngs/alea/feature/mean.hpp>
#include <alps/ngs/alea/feature/feature_traits.hpp>
#include <alps/ngs/alea/feature/generate_property.hpp>

#include <boost/cstdint.hpp>

#include <vector>
#include <cmath>

namespace alps
{
    namespace accumulator
    {
        //=================== tau proxy ===================
        //=================== tau trait ===================
        template <typename T>
        struct tau_type
        {
            typedef double type;
        };
        //=================== tau implementation ===================
        namespace detail
        {
        //set up the dependencies for the tag::autocorrelation-Implementation
            template<> 
            struct Dependencies<tag::detail::tau> 
            {
                typedef MakeList<>::type type;
            };

            template<typename base_type> 
            class AccumulatorImplementation<tag::detail::tau, base_type> : public base_type 
            {
                typedef typename base_type::value_type value_type_loc;
                typedef typename tau_type<value_type_loc>::type tau_type;
                typedef typename mean_type<value_type_loc>::type mean_type;
                typedef AccumulatorImplementation<tag::detail::tau, base_type> ThisType;
                
                public:
                    AccumulatorImplementation<tag::detail::tau, base_type>(ThisType const & arg): base_type(arg)
                    {}
                    
                    template<typename ArgumentPack>
                    AccumulatorImplementation<tag::detail::tau, base_type>(ArgumentPack const & args
                                                 , typename boost::disable_if<
                                                                              boost::is_base_of<ThisType, ArgumentPack>
                                                                            , int
                                                                            >::type = 0
                                             ): base_type(args)
                    {}
                    
                    inline tau_type const tau() const 
                    {
                        //~ //Simplebinning.h Zeile 475
                        //~ template <class T>
                        //~ inline typename SimpleBinning<T>::time_type SimpleBinning<T>::tau() const
                        //~ {
                          //~ if (count()==0)
                            //~ boost::throw_exception(NoMeasurementstag::error());
                        //~ 
                          //~ if( binning_depth() >= 2 )
                          //~ {
                            //~ count_type factor =count()-1;
                            //~ time_type er(std::abs(error()));
                            //~ er *=er*factor;
                            //~ er /= std::abs(variance());
                            //~ er -=1.;
                            //~ return 0.5*er;
                          //~ }
                          //~ else
                          //~ {
                            //~ time_type retval;
                            //~ resize_same_as(retval,sum_[0]);
                            //~ retval=inf();
                            //~ return retval;
                          //~ }
                        //~ }

                        //TODO: implement
                        return 42;
                    }
                    
                    inline void operator()(value_type_loc const & val) 
                    {
                        base_type::operator()(val);
                    }
                    inline ThisType& operator<<(value_type_loc const & val)
                    {
                        (*this)(val);
                        return (*this);
                    }
                    
                    template<typename Stream> 
                    inline void print(Stream & os) 
                    {
                        base_type::print(os);
                        os << "tag::detail::tau: " << std::endl;
                    }
                    inline void reset()
                    {
                        base_type::reset();
                    }
                private:
            };

            template<typename base_type> class ResultImplementation<tag::detail::tau, base_type> : public base_type  {

                public:

                    template<typename Accumulator> ResultImplementation(Accumulator const & accum)
                        : base_type(accum)
                    {}
// TODO: implement!
            };

        }

        //=================== call GENERATE_PROPERTY macro ===================
        GENERATE_PROPERTY(tau, tag::detail::tau)

    }
}
#endif //ALPS_NGS_ALEA_DETAIL_TAU_IMPLEMENTATION_HEADER
